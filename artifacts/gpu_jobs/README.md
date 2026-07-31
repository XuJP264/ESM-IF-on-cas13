# GPU job artifacts

Generated prediction inputs, shards, expected-output manifests, and retry
manifests live below this directory and are intentionally ignored by Git. They
can contain thousands of candidate-specific files. The tracked configuration,
export implementation, report summary, and tests reproduce them; the final GPU
bundle embeds the generated Stage-0003 input tree and verifies its hashes.
