#!/usr/bin/env python3
# Team layer over Syed Asad Ali's v12 dqi generator (Apache-2.0; see ../../NOTICE: Syed Asad Ali is the
# original builder of this decoder line). New file added by the team (Maxime Bendavid), 2026-11-18.
"""Build the v7 decoder with a configurable extractor cutoff j_min (idea DQI-A03).

Same flow as `controller_work/after_v3/integration_v7/build_v7.py`, with two differences:

1. After `Pipeline.configure()` (which runs every SHA-256 check of the pinned sources unchanged), the key-prefix
   emitter is replaced **at run time** so that the extractor receives `--j-min` instead of the hard-coded 15:
   `integrated_pipeline.emit_key_prefix` (the runtime pipeline) gets a copy of its own body with `j_min`
   threaded through, and `v7_pipeline.emit_key_prefix` forwards the team value. No pinned file is edited (the
   method Syed Asad Ali uses for his invsqrt patches). The extractor adapter in use (`first3_retained_norm_adapter`)
   accepts any 3 <= j_min < 35.
2. The three checks that pin the j_min = 15 circuit (`static_T == 384003`, equality with the j_min-15 source
   ledger, and the frozen event inventory behind the exact expected-Toffoli snapshot) become recorded
   measurements: the circuit is judged by the official verifier and `tools/gate.py`, not by the old ledger.

Usage (CPython 3.13, from the build work tree):
    python team/build_team.py --j-min 18 --output OUT.kmb --log-directory LOGS
"""
from pathlib import Path
import argparse
import gc
import json
import sys
import time

HERE = Path(__file__).resolve().parent
V7 = HERE.parent / 'controller_work' / 'after_v3' / 'integration_v7'
sys.path.insert(0, str(V7))

from v7_binding import atomic_json, digest, verify_preservation, write_inventory  # noqa: E402
import v7_pipeline as Pipeline  # noqa: E402

ORIGINAL_J_MIN = 15


def install_j_min(j_min: int) -> None:
    """Replace the key-prefix emitters so the extractor gets `j_min`. Call after Pipeline.configure()."""
    assert 3 <= j_min < 35
    P = Pipeline.P
    assert P is not None and P.__name__ == 'integrated_pipeline'

    def emit_key_prefix(builder, syndromes, *, j_min=ORIGINAL_J_MIN, stage=lambda _: None):
        # Copy of integrated_pipeline.emit_key_prefix (integration_v3) with the cutoff passed through.
        j_min = TEAM['j_min']
        km, Public, Extractor, Tail, Seed, Synthesis = P.km, P.Public, P.Extractor, P.Tail, P.Seed, P.Synthesis
        assert P._options is not None and len(syndromes) == 70
        Public.emit(builder, syndromes)
        heads = [km.array([*syndromes[2*i], *syndromes[2*i+1]]) for i in range(35)]
        carry = builder.alloc_qubits(24)
        builder.x(carry[0])
        extraction = Extractor.emit_prefix(builder, heads, carry, j_min=j_min, stage=stage)
        extraction.tail.release_spare()
        gauge = Tail.gauge_prefix(builder, carry)
        seed = Seed.emit_prefix(builder, carry[:12], heads[34])
        extraction.tail.seed_prefix(carry[:12], heads[34][12:])
        N, D = [carry[:12]], [heads[34][:12], heads[34][12:]]
        synthesis = []
        for i in range(33, -1, -1):
            Pi = [heads[i][:12], heads[i][12:]]
            rec = Synthesis.emit_prefix(builder, N, D, Pi, i, guarded=False)
            synthesis.append((N, D, Pi, rec))
            N, D = D, N + Pi
            stage(dict(stage='synthesis_prefix', node=i, degree=len(D) - 1))
        assert len(N) == 35 and len(D) == 36
        return dict(builder_id=id(builder), syndromes=syndromes, heads=heads,
                    carry=carry, extraction=extraction, gauge=gauge, seed=seed,
                    synthesis=synthesis, N=N, D=D, options=P._options, closed=False)

    TEAM['j_min'] = j_min
    P.emit_key_prefix = emit_key_prefix


TEAM = {'j_min': ORIGINAL_J_MIN}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--j-min', type=int, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--log-directory', type=Path, required=True)
    args = parser.parse_args()
    if not __debug__ or sys.version_info[:2] != (3, 13):
        parser.error('Use ordinary Python 3.13 without optimization.')
    args.output, args.log_directory = args.output.resolve(), args.log_directory.resolve()
    assert not args.output.exists() and not (args.log_directory / 'build.json').exists()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.log_directory.mkdir(parents=True, exist_ok=True)
    Pipeline.configure()
    before = verify_preservation(Pipeline.Runtime)
    install_j_min(args.j_min)
    import kickmix as km
    from v7_events import Events
    events = Events()
    events.install(Pipeline)
    started = time.monotonic()
    builder = km.CircuitBuilder()
    builder.set_max_qubits(889)
    S = [builder.create_quantum_register(12, f'S{i}') for i in range(70)]
    E = [builder.create_classical_register(12, f'E{i}') for i in range(4095)]
    stages = []

    def stage(message):
        row = dict(stage=message, allocated_Q=builder.num_allocated_qubits,
                   elapsed_seconds=time.monotonic() - started)
        stages.append(row)

    # v7_pipeline.emit keeps its own `j_min == 15` check on its argument; the team value reaches the extractor
    # through the replaced integrated_pipeline.emit_key_prefix (v7_pipeline.emit_key_prefix forwards to it).
    contract = Pipeline.emit(builder, S, E, j_min=ORIGINAL_J_MIN, scan_positions=4095, stage=stage)
    contract['j_min'] = args.j_min
    assert builder.num_allocated_qubits == 840
    circuit = builder.finish_circuit()
    del builder
    gc.collect()
    static_T = circuit.operation_statistics().num_toffoli_like
    ledger_j15 = Pipeline.source_ledger(contract['scan']['static_T'])['total_static_T']
    assert circuit.num_qubits <= 889 and dict(circuit.operation_statistics()).get('Z_POW', 0) == 0
    try:  # the event inventory is frozen for j_min = 15 (e.g. 21 'tail_or_gauge_fresh_zero_phase' events)
        expectation = events.snapshot(static_T, contract['scan']['point_tests'])
    except AssertionError as error:
        expectation = dict(skipped='event inventory differs from the frozen j_min-15 inventory',
                           detail=str(error)[:2000])
    instructions, bits, qubits = len(circuit), circuit.num_bits, circuit.num_qubits
    temporary = args.output.with_suffix(args.output.suffix + '.tmp')
    assert not temporary.exists()
    circuit.to_file(temporary, format='kmb')
    del circuit
    gc.collect()
    temporary.replace(args.output)
    Pipeline.assert_bindings()
    after = verify_preservation(Pipeline.Runtime)
    inventory = write_inventory(args.log_directory / 'source-inventory.json', Pipeline, before, after)
    report = dict(status='BUILT_UNVERIFIED_TEAM_V7_CANDIDATE', j_min=args.j_min,
                  circuit=str(args.output), sha256=digest(args.output), bytes=args.output.stat().st_size,
                  Q=qubits, static_T=static_T, static_T_j15_ledger=ledger_j15, static_T_delta_vs_j15_ledger=static_T - ledger_j15,
                  instructions=instructions, classical_bits=bits,
                  uniform_mask_source_expectation=expectation, source_inventory=inventory,
                  elapsed_seconds=time.monotonic() - started,
                  limitations=['A build is not a correctness PASS: gate with tools/gate.py (multi-salt) and the official verifier.'])
    atomic_json(args.log_directory / 'build.json', report)
    print(json.dumps({k: v for k, v in report.items() if k not in ('source_inventory', 'uniform_mask_source_expectation')}), flush=True)


if __name__ == '__main__':
    main()
