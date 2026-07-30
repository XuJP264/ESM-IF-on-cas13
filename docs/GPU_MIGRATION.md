# GPU and HPC Migration

GPU migration is manifest-based. The bundle records the git commit, environment
locks, containers, third-party commits, checkpoint/data hashes, missing assets,
configs, input shards, and expected output schemas. Large data and weights are
transferred separately with checksum-aware `rsync` or `rclone`; they are not
embedded in one giant archive.

After transfer:

```bash
bash scripts/bootstrap_gpu_node.sh
bash scripts/sync_assets.sh /source/asset-root /target/asset-root
bash scripts/verify_gpu_bundle.sh
bash scripts/launch_gpu_tmux.sh configs/gpu_node.example.yaml benchmark-experimental
```

The tmux launcher creates a unique session and run directory, logs output,
persists an exit code and pointer, and retains failed state without unbounded
retry.

