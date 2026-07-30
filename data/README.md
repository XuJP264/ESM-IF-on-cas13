# Data layout

- `raw/`: immutable upstream downloads (ignored).
- `interim/`: resumable normalized or shard-level products (ignored).
- `processed/`: analysis-ready Parquet/FASTA/MSA products (ignored).
- `manifests/`: versioned sources, licenses, schemas, hashes, and funnels.
- `experimental_structures/`: local RCSB PDB/mmCIF files (ignored).
- `fixtures/`: tiny synthetic or redistributable test inputs (tracked).

Raw data are never committed. Processors record failures and ambiguous
Cas13/direct-repeat pairs instead of silently dropping or accepting them.

