# CRC Workflow Entry

Main runnable entry point:

- `run_crc_pipeline.sh`

The script activates a local virtual environment and then runs the numbered analysis stages in order:

1. prepare QC
2. integrate and cluster
3. annotate major cell types
4. subtype CD8 cells
5. build state-aware UMAPs and bubble plots
6. refine subgroups and downstream interaction-ready summaries

This repository keeps the workflow shell plus source code, but does not include the required input data or generated result files.
