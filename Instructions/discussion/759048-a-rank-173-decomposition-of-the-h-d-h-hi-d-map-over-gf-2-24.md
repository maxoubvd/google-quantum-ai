# A rank-173 decomposition of the (H·D, H_hi·D²) map over GF(2^24)

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/759048 — copied 2026-11-15 (raw Markdown, unedited).

**Chih-Kai Wang** · 2026-11-13 14:32 UTC · 0 votes

This is a small algebraic result from my work on the DQI circuit: a 173-product decomposition of a bilinear map that appears as a repeated block in one of the public ~392k-Toffoli DQI circuits. The circuit I measured spends 176 Toffolis per block. The newest high-scoring circuits use a different structure (section 5), so I am sharing this as reusable math, with a short checker.

## 1. The map

Let K = GF(2)[t]/(t^12 + t^6 + t^4 + t + 1) = GF(2^12). Let γ be the smallest nonzero element (bit mask of t) with Tr_{K/F2}(γ) = 1 (it is 2048, i.e. t^11), and

L = K[z]/(z² + z + γ) = GF(2^24).

Write an element as a0 + a1·z with a0, a1 in K, encoded as the integer a0 | (a1 << 12). For H in L let H_hi = a1(H), which equals H + H^(2^12) and lies in K. The map is

J_0(H, D) = ( H·D , H_hi·D² ),  H, D in L, output 48 bits.

It is F2-bilinear. A decomposition J_0(H, D) = XOR_i λ_i(H)·μ_i(D)·w_i, with F2-linear forms λ_i, μ_i and 48-bit vectors w_i, runs as r Toffolis (accumulate the parity of λ_i·H and of μ_i·D on one line each with CNOTs, one Toffoli, fan out w_i). So a rank-r decomposition gives an r-Toffoli block. (This bound goes one way: blocks whose Toffoli controls mix H and D bits sit outside this bilinear form and may use fewer.)

In the circuit I examined, 182 blocks are J_0 up to linear changes of basis on H, D and the output (plus a linear map on the output lines). I counted 176 Toffolis in each. The origin of the 176 is unclear to me. For comparison, the plain interpolation route on P¹(F8), with H·D at 99 products, the two K × K products at 84 and one shared pair, gives 99 + 84 − 5 = 178. Both routes pay 99 for H·D; what the elliptic curve adds is a second aligned pair.

## 2. Construction

**First component, H·D.** L has degree 8 over F8 = GF(8). Take the curve

E: y² + xy = x³ + x² + 1 over F2.

It has 2, 14, 56, 4144 points over GF(2), GF(8), GF(64), GF(4096). So there are 14 F8-rational points (including O) and 21 places of degree 2.

Pick T in E(L) of degree 8 over F8 and a rational point R′. Use L(8·O), with basis 1, x, y, x², xy, x³, x²y, x⁴, for H, and the translate L(8·R′) (functions f(P − R′)) for D. Evaluation at T is F8-linear from each space to L, and it is a bijection when the sum A of the 8 conjugates of T satisfies A ≠ O (for H) and A ≠ 8·R′ (for D); I only use such T. Then H and D correspond to functions f_H, f_D. Their product lies in L(8·O + 8·R′), of dimension 16, and is determined by its values at the 14 rational points plus one degree-2 place (counted twice), using the leading coefficient at the poles O and R′. The degree-2 place has to be chosen so that the 16 conditions are independent (19 of the 21 places worked in my runs). H·D = (f_H f_D)(T) is then a fixed F8-linear combination of those values.

Cost: 14 F8-multiplications at 6 products plus one GF(64)-multiplication at 15 = 99.

**Second component, H_hi·D².** H_hi lies in K and D² has two K-coordinates, so this is two K × K multiplications. K has degree 4 over F8, so each takes 7 evaluation points on P¹(F8) at 6 products: 2 × 42 = 84. Computed separately, the total is 99 + 84 = 183.

**Pairs.** Since squaring is F2-linear, each second-component block has the form α(H)·β(D)² with F8-valued F8-linear forms α, β. If a first-component block uses the same forms up to F8* scalars, the two products (α·β, α·β²) can be computed as one pair block. By exhaustive search this pair has F2-rank 7 over GF(8), against 6 + 6 = 12 separately, so each aligned pair saves 5.

**Alignment.** A first-component block at point P uses α_P(H) = Tr_{L/F8}(η_P·H), where η_P is the dual element of evaluation at P. Rescaling the evaluation functionals by a on the H side and b on the D side is free (the output basis absorbs it), and this normalises one point P1 to match a second-component block exactly. A second point P2 then matches iff

η_{P2}/η_{P1} ∈ K  and  η′_{P2}/η′_{P1} ∈ √(1+z)·K

(η on the H side, η′ on the D side, the twist coming from the squaring), which has probability about 2^−12 each. I used different divisors 8·O and 8·R′ (R′ ≠ O) on the two sides so that the two dual families differ. With the same divisor on both sides η′_P = η_P up to F8* scalars, and then the two conditions contradict each other, because K and √(1+z)·K meet only in 0.

**Search.** I scanned all 15,580,848 degree-8 points T with A ≠ O (one pass over the 2^24 values of x takes about 90 s on one core). There were 1392 hits (T, R′, P1, P2). They involve 384 points T on 192 degree-8 places; T and T^(2^12) give the same hits with P1 and P2 exchanged, so there are 696 distinct configurations up to conjugation. All 1392 rebuilt to a 173-term decomposition that passed the checks below. Two aligned pairs give 99 + 84 − 2·5 = 173.

I also tried a third aligned pair at each of the other 12 rational points: the H-side condition held for 0 of 16,704 candidates (about 4 expected if it behaved randomly), and the D-side condition would cut that further. So for this curve, these divisors and alignment at rational points, 173 is the best I found. Pairing at the degree-2 place and other curves remain open.

## 3. Verification

- The map is bilinear, so equality on all 24 × 24 = 576 basis pairs and all 48 output bits is a complete check. The w_i come from solving an F2 linear system once the (λ_i, μ_i) are chosen, so this tests that the 173 products span the target.
- In addition, 200,000 random (H, D) pairs matched, and 20,000 of them were checked against a second reference formula for the multiplication.
- The checker in the appendix has its own field arithmetic, reads only the coefficient file, and repeats the 576 × 48 comparison and 20,000 random pairs. Flipping one bit of one w_i, or dropping one term, makes it fail.

## 4. Reproducing

Save the appendix files as `rank173_check.py` (numpy only, about 1 s) and `rank173_coeffs.json`:

```
$ python rank173_check.py rank173_coeffs.json
terms (= products = Toffoli): 173
basis pairs checked: 576 x 48 output bits, mismatches: 0
random (H,D) pairs: 20000, mismatches: 0
PASS
```

Format: `{"terms": [[lam, mu, w], ...]}`. lam, mu are 24-bit integers (bit i is the coefficient of bit i of H or D), w is 48-bit (low 24 bits: coefficient of H·D, high 24 bits: H_hi·D²). Term i adds w_i when parity(lam_i & H) = parity(mu_i & D) = 1.

## 5. Use and limits

- In a circuit with 176 products in each of 182 such blocks I replaced every block by this 173-term version: −3 per block, −546 Toffolis in total (static Toffoli-like count, about 0.14% of 392,344). Each replaced block was compared with the original as a map on all basis states of weight at most 2 on its 96 lines (the blocks are quadratic, so this determines them). A bit-parallel simulation of old against new circuit on 8,192 shots agreed on all lines, classical bits and phases.
- The newest public circuit family I looked at has a different structure: its 1,376 GF(2^24) multiplications each act on their own pair of operands, while a pair block needs two products of the same operands. So the construction applies to the older family.
- For the other block type in the older circuits, (H·D, (cH + H_hi)·D²) with c ≠ 0, one aligned pair gives 99 + 99 − 5 = 193 (checked on one example), above the 186 I know for that block. Three aligned pairs would give 183; that search (roughly 10^10–10^11 trials, each a few hundred field multiplications) is still open.
- Whether 173 is optimal is open; the alignment search covered this one curve and these divisors.

The ingredients are standard. Interpolation on curves goes back to Chudnovsky–Chudnovsky; the asymmetric version (different divisors on the two sides) is Randriambololona's (arXiv:1107.0336), and elliptic-curve algorithms of this type are worked out in Ballet–Bonnecaze–Tukumuli (arXiv:1303.7082). For small bilinear maps such as the GF(8) and GF(64) multiplications used here, see Barbulescu–Detrey–Estibals–Zimmermann (ePrint 2012/110). As far as I know, the new part is aligning the evaluation blocks of H·D with those of H_hi·D² so that two of them become rank-7 pairs. If it already appears somewhere, please point me to it.

## Appendix A: rank173_check.py

```python
#!/usr/bin/env python3
"""Check that 173 rank-one terms reconstruct J_0(H,D) = (H*D, H_hi*D^2) over L = GF(2^24). Needs numpy only.
Usage: python rank173_check.py [rank173_coeffs.json]"""
import sys, json, random
import numpy as np

KPOLY = 0x1053                      # t^12 + t^6 + t^4 + t + 1 (K = GF(2^12))

def kmul(a, b):                     # carry-less multiply, reduce mod KPOLY
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
        if a & 0x1000:
            a ^= KPOLY
    return r

def ktrace(a):                      # absolute trace K -> F2
    s, x = 0, a
    for _ in range(12):
        s ^= x
        x = kmul(x, x)
    return s

# sanity: KPOLY is irreducible (t has order 4095, so K really is a field)
x, order = 1, 0
while True:
    x = kmul(x, 2); order += 1
    if x == 1: break
assert order == 4095, "KPOLY is not primitive"
GAMMA = next(g for g in range(1, 4096) if ktrace(g) == 1)

def lmul(x, y):                     # L = K[z]/(z^2 + z + GAMMA); element = a0 | a1 << 12
    a0, a1, b0, b1 = x & 0xFFF, x >> 12, y & 0xFFF, y >> 12
    m11 = kmul(a1, b1)
    return (kmul(a0, b0) ^ kmul(m11, GAMMA)) | ((kmul(a0, b1) ^ kmul(a1, b0) ^ m11) << 12)

def J0(H, D):                       # 48-bit int: low 24 = H*D, high 24 = H_hi * D^2
    return lmul(H, D) | (lmul(H >> 12, lmul(D, D)) << 24)   # H_hi = a1(H), an element of K inside L

def bits(v, n):
    return np.array([(v >> i) & 1 for i in range(n)], dtype=np.int64)

def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "rank173_coeffs.json"
    terms = json.load(open(path))["terms"]
    r = len(terms)
    lam = np.array([bits(t[0], 24) for t in terms])         # r x 24
    mu = np.array([bits(t[1], 24) for t in terms])          # r x 24
    w = np.array([bits(t[2], 48) for t in terms])           # r x 48
    assert all(t[0] and t[1] for t in terms), "zero form"
    # bilinear => it suffices to match all 24*24 basis pairs (e_a, e_b) in all 48 output bits
    recon = np.einsum('ia,ib,ik->abk', lam, mu, w) % 2      # 24 x 24 x 48
    target = np.array([[bits(J0(1 << a, 1 << b), 48) for b in range(24)] for a in range(24)])
    bad = int((recon != target).sum())
    print(f"terms (= products = Toffoli): {r}")
    print(f"basis pairs checked: 576 x 48 output bits, mismatches: {bad}")
    # independent random test straight from the definition (no tensor, parity of masks)
    rng, nbad, n = random.Random(12345), 0, 20000
    for _ in range(n):
        H, D = rng.getrandbits(24), rng.getrandbits(24)
        o = 0
        for l, m, ww in terms:
            if bin(l & H).count("1") & 1 and bin(m & D).count("1") & 1:
                o ^= ww
        nbad += o != J0(H, D)
    print(f"random (H,D) pairs: {n}, mismatches: {nbad}")
    ok = bad == 0 and nbad == 0 and r == 173
    print("PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)

if __name__ == "__main__":
    main()
```

## Appendix B: rank173_coeffs.json

```json
{"description":"Rank-173 decomposition of J_0(H,D)=(H*D, H_hi*D^2) over L=GF(2^24)=K[z]/(z^2+z+gamma), K=GF(2)[t]/(t^12+t^6+t^4+t+1), gamma=smallest g in 1..4095 (as bit mask of t) with Tr_{K/F2}(g)=1. Element a0+a1*z is the int a0|(a1<<12). Each term [lam,mu,w]: 24-bit masks lam,mu and 48-bit output mask w; J_0(H,D)=XOR_i parity(lam_i&H)*parity(mu_i&D)*w_i, output low 24 bits = H*D, high 24 bits = H_hi*D^2.","K_poly":4179,"rank":173,"terms":[[12117590,11272192,7448380],[2477129,1294336,10942246],[8390127,4411392,13172780],[10865062,5279744,14110746],[10299935,12566528,12101392],[3728313,15683584,7337738],[1975678,9806004,3766127],[4725057,8743482,13795545],[15538563,1009049,1501284],[10813634,9047971,15464374],[5651519,1100430,3120907],[15940861,10143021,12871869],[15138186,4179439,14742299],[11023989,5855654,3845164],[9903430,978928,14617132],[4139827,5748310,14311223],[5164031,6724681,4192567],[7463116,3222047,15050240],[13461476,3792362,5897531],[2260320,12102866,15780980],[1396602,10818653,8010264],[3617306,1947791,11088207],[15669892,8483128,2344739],[14166174,10275255,9107052],[3519979,6145618,12132762],[12735030,16619481,12543983],[7292459,6188084,3208180],[11342877,10746861,411253],[16246749,10506635,9032302],[5960640,240230,9409563],[12116534,15888779,6282092],[2474944,16498111,5863964],[8389661,8857561,8266656],[10864605,8165990,434544],[10297846,642100,2226380],[3728939,7689810,2578876],[13307524,7419906,14292295],[15565008,1542036,15684527],[5739913,13780287,1317808],[12195161,12960427,3491560],[2526804,6733718,13504247],[10263309,10714429,16467999],[5748310,11779964,12910499],[6724681,13493314,17273],[4179439,14932168,3078420],[5855654,3030154,12893402],[3222047,8280894,15337143],[6845369,5269428,3062381],[1912442,6405730,13587417],[12119889,11608714,6537401],[3729348,13275321,15454489],[8390805,8106547,11309920],[10862891,13671656,2392768],[2476478,11218651,8917408],[15742396,12656547,13340208],[1959294,1247886,16079405],[8637693,15028404,3378072],[10365315,16144954,4117533],[15585474,13767981,16254376],[7600449,2379543,13029813],[10639300,3371193,5842790],[10578069,10695219,15826530],[6604527,5782609,4393244],[12955258,16453218,11032836],[212817,9454218,1715834],[13013291,7032040,11695998],[8251364,12904887,6133862],[8940896,1414237,16314349],[5767034,7553336,2107022],[14651930,6739301,10844043],[16091780,13729258,8240872],[2758814,12036239,14207331],[12260599,2558225,3152563],[12022389,5146741,16468920],[812674,6914404,13324555],[7425434,8969195,2876351],[3201218,13950195,1132933],[4298072,6031128,3844666],[15806689,4172190,12893757],[8375835,12461624,15668617],[9371386,8490918,2862004],[340687,7494212,333748],[513996,2959091,14307461],[189699,6257847,14632753],[8736692,12920881,352191],[4768021,4709432,1025853],[13471393,8582153,718978],[9015296,16273408,267524922933248],[14262272,13844480,40338332844032],[4284416,6971392,95657511616512],[10010624,12132352,237357072646144],[5255168,2822144,180663504338944],[13168640,9580544,127199751438336],[3272704,11366400,207876417126400],[9334784,7647232,114074331381760],[16388096,4132864,16149077032960],[7626752,4956160,240586888052736],[12550144,14270464,197499776139264],[13361152,9592832,115929757253632],[3112960,4243456,134690174402560],[3465216,11640832,22265110462464],[11595776,14233600,89335319756800],[8654848,6852608,121770912776192],[1794048,15818752,48103633715200],[10448896,10088448,75866302316544],[3235840,90112,61641370632192],[12050432,5627904,258179074097152],[16457728,13320192,25632364822528],[5029888,10395648,231722075553792],[8814592,5537792,51951924412416],[13254656,13246464,278726197641216],[704512,11542528,280512904036352],[8261632,3379200,67482526154752],[3031040,11837440,140737488355328],[5263360,8859648,213580133695488],[7655424,8630272,139775415681024],[2392064,294912,208220014510080],[11505664,10801152,162384123527168],[14327808,11993088,169256071200768],[9932800,14536704,141424683122688],[5058560,7000064,11269994184704],[7680000,1298432,20959440404480],[3670016,7929856,27831388078080],[8388608,8388608,168156572955873],[3727360,10297344,26182125284864],[2473984,12115968,198186985340382],[12115968,3727360,49134434966241],[1908736,10862592,65283513045982],[10297344,1908736,157367607420735],[10862592,2473984,179838871669055],[10424320,7013387,33403437056],[13987840,15517544,50532974592],[13107200,9060212,29091692544],[1929216,6748188,51589939200],[4874240,8897379,4311744512],[5705728,14760831,55868129280],[11132928,7650167,53166997504],[6213632,9595874,1577058304],[15663104,4957258,49476009984],[11653120,14274472,52596572160],[16199680,15126677,33889976320],[4644864,4136765,50918850560],[720896,4973874,40869298176],[1744896,16611779,11525947392],[12853248,1326992,9076473856],[14581760,15287891,47966060544],[1155072,11967729,49811554304],[13574144,6281890,2986344448],[16474112,908459,15988686848],[16105472,16342846,57747177472],[14684160,4036312,34426847232],[1429504,12896742,63535316992],[958464,16024469,50415534080],[1798144,3164787,23454547968],[6545408,11497377,65716355072],[1896448,12808991,38218498048],[16003072,2742503,8975810560],[15253504,15379448,32933675008],[8327168,7085246,57344524288],[9949184,8828742,46992982016],[16093184,9561997,38721814528],[15331328,5468716,402653184],[11673600,5809107,61773709312],[6017024,774655,39124467712],[1859584,12752289,31776047104],[4698112,13190238,62176362496],[13885440,12293049,35654331896],[10149888,477606,47761427250],[1511424,4714480,11798986417],[4784128,16546889,47024727754],[9302016,16463343,8263546243],[6164480,12376607,14288050299],[12906496,4236886,41363166025]]}
```
