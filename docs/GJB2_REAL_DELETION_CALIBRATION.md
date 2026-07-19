# GJB2 real-deletion validator calibration

A bounded, receipt-bearing calibration of the MCORE-1 metrical-tree validator on
two **real** GJB2 single-base deletions, plus a committed-WAV codec round trip.

## Experiment question

> Do the two real GJB2 deletions `c.35delG` and `c.235delC` produce
> deterministic, inspectable validator error geometries under the current
> MCORE-1 metrical-tree checker, and does that geometry survive the repository's
> Gabor WAV codec round trip?

It does **not** ask whether MCORE-1 explains pathogenicity, whether an encoder
carry causes the divergence, or whether LLM hallucinations are deterministic.

## The two real variants

Both are common, real GJB2 (Connexin 26) coding-sequence deletions modeled
against RefSeq **NM_004004.6** (681-bp CDS):

- **`c.35delG`** — deletes `G` at CDS position 35 (1-based).
- **`c.235delC`** — deletes `C` at CDS position 235 (1-based).

Each yields a 680-nt mutant. This is a **real-data calibration**, not broad
external validation: two positive single-base deletions in one gene.

## Repository boundary

| Concern | Repository |
|---|---|
| DNA→trit encoder, metrical-tree validator (`check_deletion`), generic deletion-shape schema + signatures, generic tests | `mcore-1` |
| NM_004004.6 reference handling, `c.35delG`/`c.235delC` construction, Gabor WAV codec, calibration runner, receipts | `gjb2-mcore-sonification` |

`mcore-1` is imported as the **source of record** (never copied in); the loaded
`mcore_1.__file__` is recorded and validated so a stale installed package cannot
be used silently.

## Zero-carry invariant (why this is not a carry effect)

The DNA mapping is `A=0, C=1, G=2, T=0` with a T-bias `εT=1`, and per position
`u = v + carry + εT`, `t = u mod 3`, `carry = ⌊u/3⌋`.

From `carry = 0`, every base yields `u ≤ 2`, so `carry_out = 0`. Zero is an
absorbing fixed point: **no string over `{A,C,G,T}` ever produces a nonzero
carry.** The encoder is therefore *state-compatible but carry-inert on valid DNA
input*. This is proved as an executable invariant in
`mcore-1/tests/test_carry_invariant.py`. Consequently any observed deletion
error geometry **must not** be described as a DNA-encoder carry effect. The
carry machinery is retained; only the present invariant is made explicit.

## Deletion-shape schema (`mcore_1.deletion_shape`)

`summarize_deletion_shape(rows, *, deletion_pos_1, wt_length, mutant_length)`
turns `check_deletion` rows into a canonical, hash-stable `DeletionShape`.

Key properties:

- **UUID-independence.** `NodeResult.node_id` is a random `uuid.uuid4()` that
  changes every run; it is **excluded** from all canonical output. Node identity
  is the deterministic span `(leaf_lo, leaf_hi)` plus post-order rank, so the
  canonical JSON and both signatures are byte-identical across repeated calls.
- **Three signatures.**
  - `receipt_signature` — SHA-256 over the **complete absolute** artifact
    (absolute spans, `k`, counts): the identity of a specific deletion's receipt.
  - `topology_signature` — SHA-256 over the geometry **relative to `k`** (offsets
    `lo−k`, `hi−k`) but **retaining each node's post-order rank** in the whole
    mutant tree. Because that rank counts nodes to the left of the deletion, this
    digest is position-dependent through placement — use it for exact topological
    identity.
  - `geometry_signature` — SHA-256 over a **truly local, position-stripped** view:
    a canonically sorted multiset of `k`-relative node descriptors (offsets,
    widths, survivor counts, site class, error kinds) with **no post-order rank
    and no post-order ordering**. Two deletions with the same local error geometry
    share this digest regardless of where they sit in the sequence.
  - Local-geometry equality between deletions is tested with `geometry_signature`
    **only** — never the receipt signature; use `topology_signature` when exact
    tree placement must also match.
- **Coordinate semantics.** `coordinate_width = hi − lo + 1` (original-coordinate
  hull, may include the deleted site) versus `survivor_leaf_count` (surviving
  original indices under the node, never counting `k`); they differ exactly when
  a hull straddles `k`.
- **Determinism / robustness.** Output depends only on the *set* of node records
  plus `(k, wt_length, mutant_length)` — input row order is normalized
  deterministically via a reconstructed post-order. Duplicate spans, spans not
  matching the topology, `valid`/`errors` inconsistencies, and degenerate spans
  raise `ValueError`. All-valid trees yield `null` first/narrowest/widest spans.

Signed scientific content excludes timestamps, machine paths, Git SHAs, and
environment metadata; the schema version is included in the signed payload.

## Prefix-index vs gap-restored controls

Two minimal linear controls separate the validator's tree geometry from a raw
displacement artifact:

- **Prefix-index (deletion-induced suffix shear, displacement-sensitive).**
  Compare WT and mutant trits at the same array index after the deletion. Every
  position at/after the cut is shifted, so this diverges wherever consecutive WT
  trits differ.
- **Gap-restored (deletion-localized).** Insert a gap sentinel at the deleted WT
  column and compare surviving positions only. Because the encoder is
  carry-inert, **surviving trits are unchanged → zero divergence.**

The gap-restored zero therefore shows the validator's CONSERVATION/OVERFLOW
outcomes come from **frozen-topology re-pooling**, not per-position trit edits;
the prefix-index count is the shear artifact and must not be read as the tree
geometry.

## WAV round-trip method

The committed `audio/` files encode one 40 ms Gabor frame per trit at 48 kHz,
mono, 16-bit, single carrier per trit (`0→800, 1→1600, 2→3200` Hz). The
pure-stdlib decoder (`code/wav_decode.py`) validates format, exact frame
divisibility, and expected lengths, then recovers each trit by Goertzel energy
argmax with a defined per-frame confidence = normalized winner margin over
spectral energies `(E_win − E_runner_up)/E_win ∈ [0,1]` (these are Goertzel
powers, not probabilities).

The committed artifacts are the WT stream (681 frames) and per-allele **delta**
streams (680 frames). The round trip decodes WT and delta **independently**,
checks the decoded delta against the WT-derived delta, and reconstructs the
mutant `mut[i] = (wt[i] + delta[i]) mod 3`, which must equal WT with the deleted
column removed. Reconstructed WT/mutant streams then feed `check_deletion`, and
the resulting deletion-shape signatures are compared to the DNA lane when a
reference is present. Source WAVs are never overwritten.

## Results

Run:

```bash
python code/real_deletion_calibration.py            # audio lane + fail-closed DNA gate
python code/real_deletion_calibration.py --fasta data/refseq/NM_004004.6.fasta   # full run
```

Receipts are written under `artifacts/calibration/gjb2-real-deletions/`
(`calibration.json`, `deletion_shapes.jsonl`, `node_results.csv`, `summary.md`,
`run_manifest.json`, `implementation_*.patch`). The manifest hashes every other
artifact (and never itself), records repo branches/HEADs, dirty status, Python
and dependency versions, reference provenance (raw + normalized CDS hashes;
`hash_recorded` vs `hash_verified`), the exact commands, and per-criterion
pass/fail.

Environment note: NCBI E-utilities egress is blocked where this runs, so the DNA
lane is **fail-closed** — it halts at the reference gate unless a verified FASTA
is supplied (`--fasta` / `$GJB2_CDS_FASTA` / `data/refseq/NM_004004.6.fasta`).
The audio lane runs regardless and is labeled **"committed audio-artifact-derived;
biological-reference cross-check pending."**

Enforced invariants (each forces a non-zero exit and is recorded in the manifest,
distinct from a soft reference-gate halt): a *supplied* reference failing its
guards or a supplied `--expected-cds-sha256` mismatch; unstable signatures across
repeats; upstream/local encoder disagreement or nonzero carry; WAV codec
self-inconsistency; and, when the DNA lane runs, DNA↔audio signature mismatch. A
missing reference is not a failure — it is the expected fail-closed halt.

On the committed audio artifacts (WT decode margin 1.0000, codec round-trip
self-consistent):

| allele | k | invalid nodes | geometries stable (3×) |
|---|---|---|---|
| `c.35delG` | 35 | 595 | yes |
| `c.235delC` | 235 | 599 | yes |

The two alleles' `geometry_signature`s **differ** (distinct error geometries).
Consult the committed `run_manifest.json` for the exact signature digests of the
current run.

## Claim tiers

- **Established in code:** deterministic DNA→trit mapping; strict
  681/`ATG`/`G@c.35`/`C@c.235` reference guards with no silent fallback;
  validator `NodeResult` output; carry-inert DNA mapping (executable invariant);
  repeatable deletion-shape receipt + geometry signatures.
- **Supported by this calibration:** whether the two real deletions produce
  stable, inspectable error geometries, and whether the committed Gabor WAV codec
  preserves them through decode and revalidation.
- **Not established:** biological mechanism; clinical utility / pathogenicity
  prediction; broad external generalization; LLM-hallucination determinism.

## Limits

- Two positive single-base deletions in one gene — a calibration, not validation.
- No clinical or mechanistic claim; the encoder is a deterministic lexicon.
- The observed geometry is a frozen-topology re-pooling effect, **not** a carry
  effect.
- The audio lane is derived from committed artifacts; the biological-reference
  cross-check is pending a supplied NM_004004.6 FASTA.
