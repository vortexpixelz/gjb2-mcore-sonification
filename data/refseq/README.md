# RefSeq reference drop-in for the real-deletion calibration

The real-deletion calibration (`code/real_deletion_calibration.py`) is
**fail-closed**: its DNA lane will not run against the embedded demonstration
fallback. It requires the verified NM_004004.6 reference.

Drop the reference here as:

```
data/refseq/NM_004004.6.fasta
```

Accepted forms (the loader distinguishes and records which):

- the **full NM_004004.6 transcript** FASTA — the loader extracts the CDS as
  `seq[178:859]` (GenBank CDS `179..859`), or
- the **681-bp CDS** directly (must start `ATG`).

The loader verifies length 681, `ATG` start, `G` at c.35, and `C` at c.235,
records both the raw-file SHA-256 and the normalized-CDS SHA-256, and marks
`hash_verified` only if you pass an expected hash
(`--expected-cds-sha256` / `$GJB2_CDS_SHA256`).

You can also point elsewhere without copying a file here:

```bash
python code/real_deletion_calibration.py --fasta /path/to/NM_004004.6.fasta
# or
GJB2_CDS_FASTA=/path/to/NM_004004.6.fasta python code/real_deletion_calibration.py
```

NCBI E-utilities egress is blocked in the calibration's execution environment,
which is why the reference is supplied as a hashed file rather than fetched live.
