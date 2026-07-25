# GJB2 real-deletion calibration — committed receipt index

Compact, in-tree provenance record for the full DNA + audio calibration.
The large generated receipts (`calibration.json`, `deletion_shapes.jsonl`,
`node_results.csv`, `run_manifest.json`, `implementation_*.patch`) are
**gitignored** under `artifacts/calibration/gjb2-real-deletions/` and are
reproducible via the command below; their SHA-256s are recorded here.

## Exact command
```bash
python code/real_deletion_calibration.py \
  --fasta data/refseq/NM_004004.6.fasta \
  --expected-cds-sha256 4e200a0cd3e11879057fe0a2557e25de6925934e798e63ca4dc9235dec08907a
```
Process exit code: **0** · `failed_mandatory_criteria`: **none** · `invariant_violations`: **none**

## Reference
- accession: **NM_004004.6** · input_kind: `cds` · extraction: `identity`
- normalized CDS SHA-256: `4e200a0cd3e11879057fe0a2557e25de6925934e798e63ca4dc9235dec08907a` (**hash_verified=True**)
- raw FASTA file SHA-256: `1bdb7716f21e366d870754f6ea0a8dcae23d353c216a5f32ded805d16b2aa72c`
- guards: {'length_681': True, 'starts_atg': True, 'g_at_c35': True, 'c_at_c235': True}

## Commit SHAs (run provenance)
- gjb2:   `e34176ea3d45865dc1f6235c6eb333c160aa7bd2` (branch `claude/code-work-order-execution-1jzzr6`, dirty=False)
- mcore-1: `14252ebccf3283c42a845dc01a166e3c8cac65a6` (branch `claude/code-work-order-execution-1jzzr6`)

_This index is committed in the immediately following commit; the receipts above were generated at the gjb2 SHA shown._

## Acceptance criteria
| criterion | status | mandatory |
|---|---|---|
| `mcore1_existing_tests_pass` | verified_out_of_band |  |
| `deletion_shape_tests_pass` | verified_out_of_band |  |
| `carry_inert_invariant_test_passes` | verified_out_of_band |  |
| `prefix_and_gap_controls_emitted` | pass |  |
| `report_not_attributing_to_carry` | pass |  |
| `artifacts_hashed_in_manifest` | pass |  |
| `committed_wavs_not_overwritten` | pass |  |
| `strict_reference_guard_no_fallback` | pass | yes |
| `three_runs_identical_hashes` | pass | yes |
| `both_variants_canonical_shape` | pass |  |
| `encoder_streams_match_upstream_local` | pass | yes |
| `carry_logs_zero` | pass | yes |
| `wav_codec_self_consistent` | pass | yes |
| `dna_audio_cross_check_equal` | pass | yes |
| `dna_audio_raw_streams_equal` | pass | yes |

## Generated artifact SHA-256 (gitignored; reproducible)
| file | sha256 |
|---|---|
| `calibration.json` | `580bebbaa7d508aaf9b9bc0251970d7dea5227ead31531f35b25780ab046e10b` |
| `deletion_shapes.jsonl` | `ff22b29a2858d5fc4528818ca5ed235d89ee87ebde240fdab900756209a7f9d0` |
| `node_results.csv` | `80d4e6ba5fdccb4cb27ce26d83672c884aaf950bbdb9918f928baa86c274b6b1` |
| `summary.md` | `0ea20dcc72381bd7303e44a07971024807e28fdd508e643010eae79243ddc3b9` |
| `implementation_gjb2.patch` | `88b003146103d565742c101754c2d80c7c093fa6b1afeef30d70ffef7399d72a` |
| `implementation_mcore1.patch` | `4353941b5140223e7fe6ee73c3a43498a19bcf48162641168a54024eb870eb87` |
| `run_manifest.json` | _(not self-hashed; contains a run timestamp)_ |

## Raw trit-stream cross-check (DNA vs audio) — exact
| stream | length | equal | mismatch | first_mismatch | shared SHA-256 |
|---|---|---|---|---|---|
| `wt` | 681 | True | 0 | None | `6bf1b73341d4235adf29b8fd20d0d65bc776746bb41fb0cead29d5705a6e89a4` |
| `c35delG_delta` | 680 | True | 0 | None | `27c2974d6c2e2b778ef4d4c6140047f9987e6e9f75025bb8cd95c53f88edd3c6` |
| `c35delG_mut` | 680 | True | 0 | None | `cf693de7d04de8ac004ec52fcfe03ed8d8bdccce5f3cb36e24cfa435a69e5b44` |
| `c235delC_delta` | 680 | True | 0 | None | `9e30cf00ac87ced9354439c9377664171113b45f6d88deac9506cc0927874093` |
| `c235delC_mut` | 680 | True | 0 | None | `061a8b7e1b75c740faf0fb1670e38cf3659725d5a503ab99b3774f69fd267551` |

## Non-claims
- The observed geometry is NOT attributed to a DNA-encoder carry effect (mapping is carry-inert).
- Two variants are a real-data calibration, NOT broad external validation.
- No biological, diagnostic, or clinical claim is made.
