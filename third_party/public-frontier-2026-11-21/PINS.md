# Public frontier set of 21 Nov 2026 (board 33.10186)

Third submission candidate (milestone M2.E). None of these circuits is ours. Source: the output of goodpjw2008's public
notebook "QCC | Dead-Word Loader Chain DP | LB 33.10186" (Apache-2.0); mizeroluckygall's re-run ships the same members.
Downloaded 2026-11-21 with `kaggle kernels output` into `downloads/outputs-2026-11-21/` (git-ignored, read only).

Lineage (from `SOURCE_CREDITS.md`, verbatim copy here; Syed Asad Ali first as his NOTICE requires):

- invsqrt `d2148a7c` (37 qubits): Syed Asad Ali (generator line) → Team ThemBones (Danish Ali) → Almond Labs → Rithvik Redrouthu.
- dqi `4b572141` (900 qubits): Syed Asad Ali → Team ThemBones → Tony Li (tonylica).
- dfthcsos: table `fc4b4c01` by Ertugrul Ozer (ozertuu), refining Almond Labs' `6b9903dd`, words and lineage from Team ThemBones;
  loaders by Rithvik Redrouthu, rupam2020 (No Bodies) and boristown (team 暗黑AGI), chained by goodpjw2008's dead-word DP.

Licences: Apache-2.0 for every part (authors' statements); `SYED_ASAD_ALI_NOTICE` copied verbatim.

Zip: `submission.zip`, 77,981,775 bytes, SHA-256 `96fa5226d56d937a462777d1444e01391bc6f5c72e1c8859d1e711955d2710a2`.

Costs as reported in the author's `receipt.json` (official verifiers, default salt):

| Task | Qubits | Toffolis | Cost |
|---|---:|---:|---:|
| invsqrt | 37 | 1406.31 | 60,955.33 |
| dqi | 900 | 316131.2 | 3.050416e+11 |
| dfthcsos (brot 7, λ_eff 46.91821) | 405 | 13290 | 252,534,924.82 |

| Member | SHA-256 |
|---|---|
| `invsqrt.kmb` | `d2148a7cd7a4d0a004eca1e7b18f9fd1282df4ed3be4a9aa970842f9119fa307` |
| `dqi.kmb` | `4b5721411adc239f61fc11c773a2358b75dfe7cb15586efae4f0459c77d77e52` |
| `dfthcsos/angle_k_0.kmb` | `8df36485e020b9e30353142a694b69018dc7defdf81847a6b405cd0d9b152402` |
| `dfthcsos/angle_k_1.kmb` | `b84a398bc050db8dd8b4caa5a73ebee8a1ff262453ce9397671775718863f492` |
| `dfthcsos/angle_k_2.kmb` | `046679fe31313ac2393df03c85e83f3e9a313909919a3288c0d4f109e09951e2` |
| `dfthcsos/angle_k_3.kmb` | `9ee7708d0b30404576ed7617cadd9e842228286d54104de54b1c2fde4437a575` |
| `dfthcsos/angle_k_4.kmb` | `cf18699e46b45367ff88d8d31140543738a2e3a7dcb559c5497b9c7b53a0f9d4` |
| `dfthcsos/angle_k_5.kmb` | `0718ed8e096bb8947de07341240448ba2b1ae399801bbc9426ab5e9fa8d789b4` |
| `dfthcsos/angle_k_6.kmb` | `64c95e3cf1be5e941b4f6e39a811cf11e064f355b40d65c61ea5c7152b9c48e3` |
| `dfthcsos/angle_k_7.kmb` | `18da4cf5c0c67ca2f169abea481c1aa6ba26db9d3e2fb79263267f7ca5d65dff` |
| `dfthcsos/angle_k_8.kmb` | `0d7574240d1dbf4a6af8fcdfea41a3284dc460bb62451a0320c0a6e4ebce038a` |
| `dfthcsos/angle_k_9.kmb` | `a55398d3676b1710a48a3a23e9e98fb6e89b3732a47d17e9bb9b7b8a5de7cf12` |
| `dfthcsos/angle_k_10.kmb` | `31d55da61c233c08be359a1bb5fb61c79db5891bb840188783191a382f5c9c50` |
| `dfthcsos/angle_k_11.kmb` | `5ac10d8f786b6e54f161ee52c7b6a3b42f6b3b76b06edbb9bb248f3f84217097` |
| `dfthcsos/angle_k_12.kmb` | `af769644c2241f34fa430106e1bb357adce44e12c2a38dbf4267f80aeff258c6` |
| `dfthcsos/angle_k_13.kmb` | `fdc1561ef8347941f7db029bfb1b83ce661405b7b4137f47f1b6b6a72ea2e849` |
| `dfthcsos/angle_k_14.kmb` | `b2105608105224c01e3198c3d382c7d6294e37f7fbbecf5cd9c26602a9d19ebf` |
| `dfthcsos/submission.npz` | `fc4b4c012372517daff3e5c42062ae3ff072848f165d768d4468b6c4c4ac76b6` |
