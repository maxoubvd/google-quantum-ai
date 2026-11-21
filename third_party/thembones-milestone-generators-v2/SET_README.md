# Verified three-task portfolio (team ThemBones)

Circuits for invsqrt, DQI and dfthcsos, checked by sha256 (MANIFEST.sha256) and packed by the companion notebook.

- `set.zip` (or the unpacked files): `invsqrt.kmb`, `dqi_fileA.kmb`, `dfthcsos/angle_k_0..14.kmb`, `dfthcsos/submission.npz`.
- `dqi_fileA.kmb` is the DQI circuit WITHOUT its re-roll label; each leaderboard ticket is this file plus one INFO_TAG
  `reroll:<label>` (labels: r46zd20Aa, r46zd20Af, r46zd20Ag); the notebook rebuilds every ticket and checks its sha256.
- `costs__ticket_<label>__<salt>.json`: the official verifiers' outputs for every ticket on the local salts "" and "alt-salt-1".
- `LICENSE` (Apache-2.0), `NOTICE`, `third_party__*` (upstream LICENSE / NOTICE / CREDITS). The folder is flat on
  purpose (the Kaggle CLI skips sub-folders by default).
- The competition data is not included.
