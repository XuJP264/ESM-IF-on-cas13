# Experimental Cas13 structure data card

Generated from locally hashed RCSB PDB/mmCIF files and RCSB Data API metadata.
This is real structural QC (Evidence Level 0), not evidence that a designed
sequence is an effective Cas13.

| PDB | Subtype | State | Method | Resolution (Å) | Protein chains | RNA chains | Design chain | Coordinates/SEQRES | Mapped identity | Chain breaks | Status |
|---|---|---|---|---:|---|---|---|---:|---:|---:|---|
| 6E9F | VI-D | crRNA_target_ternary | ELECTRON MICROSCOPY | 3.3 | A | B,C | A | 864/954 | 1.0000 | 4 | included |
| 5XWP | VI-A | crRNA_target_ternary | X-RAY DIFFRACTION | 3.086 | A,B | C,D,E,F | A | 1125/1160 | 1.0000 | 8 | included |
| 6E9E | VI-D | crRNA_binary | ELECTRON MICROSCOPY | 3.4 | A | B | A | 863/954 | 1.0000 | 5 | included |
| 5XWY | VI-A | crRNA_binary | ELECTRON MICROSCOPY | 3.2 | A | B | A | 1117/1159 | 1.0000 | 3 | included |

The primary benchmark uses 6E9F and 5XWP. Their same-study binary states, 6E9E
and 5XWY, are retained for conformational comparison. RNA atoms are preserved
for contact annotation and LigandMPNN, but are never supplied as protein chains
to ESM-IF1. Coordinate gaps are recorded rather than imputed.

Related structures were queried and recorded in the manifest, but were not
silently added to the preregistered bootstrap benchmark.
