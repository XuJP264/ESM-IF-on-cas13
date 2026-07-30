# References

`manifest.yaml` is the source of truth for bibliographic metadata, access
status, licenses, local PDF paths, and checksums. `fetch_references.sh` downloads
only entries explicitly marked open access from an official repository. The
entire `papers/` directory is ignored by Git.

Paywalled works retain title, DOI, URL, BibTeX, and notes only. Do not bypass
access controls or commit PDFs.

