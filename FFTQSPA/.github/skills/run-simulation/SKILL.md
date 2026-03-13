---
name: run-simulation
description: This skill provide the workflow for runing a simulation.
---
Each simulation run consists of the following steps:
- Change the channel parameters according to the user's needs. The channel parameters are currently hardcoded in `main.py` and passed to `DNAChannel(...)`. Other Channel parameters are listed in `python/Model/config.py`.

- Everytime the parity file is changed, the mapping file should also be updated accordingly. A n ary code with q_ary = n^2 should have exactly n signals in the mapping file. The mapping file are named as Signal_Set_{n}.txt, and they are located in the `Mapper` folder.
