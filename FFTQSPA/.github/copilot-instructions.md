# Copilot Instructions for FFTQSPA Project

## Project Purpose
This repository implements a **q-ary LDPC codec** with **FFT-QSPA decoding** in C++, and exposes the core encoder/decoder to Python via `pybind11`. 
The project also includes a simulation pipeline for evaluating the codec's performance over a DNA storage channel to compare the q-ary LDPC code with the 5G NR binary LDPC code, with results output in CSV format for BER/FER analysis. 
Then results can be plotted using the #`Data2Figure.py` script.

Typical workflow in this repo:
1. Build/load `fftqspa` Python extension from C++ sources.
2. Generate or load LDPC parity/check data files.
3. Encode binary information bits (`encoder4bibo`).
4. Simulate channel outputs (AWGN or DNA channel pipeline).
5. Decode using FFT-QSPA (`decode4bibo`) with bit probabilities `P(bit=0)`.
6. Compute BER/FER and persist CSV results.
7. Plot results with `Data2Figure.py` or similar.

---

## High-Level Architecture

### 1) C++ core (algorithm)
- `BCJRQSPA.{h,cpp}`
  - Main class exposed to Python.
  - Handles:
    - Parity/matrix loading from `.dat` file
    - q-ary LDPC encoding (`encoder4BiBo`)
    - FFT-QSPA decoding (`FFTQSPA4BiBo`)
  - Uses finite-field arithmetic and Tanner graph message passing.

- `QaryLDPC.{h,cpp}`
  - `ldpc_edge` and `Tanner_Graph` data structures.
  - Check/variable message buffers (`m_v2c`, `m_c2v`) and transforms/permutations.

- `FiniteField2.{h,cpp}`
  - GF(2^m) operations (`Add`, `Mult`, `Div`, etc.).

- `Mapper.{h,cpp}` and signal-set txt files
  - Modulation mapping support.
  - **Important consistency rule**: number of signals must equal `q_ary` in the parity file.

- Supporting files: `Random.*`, `Interleaver.*`, `Qary_Gauss_E.*`, `lsfr.*`, `util.*`, `stdafx.h`.

### 2) Python bindings and simulation
- `python/fftqspa_bindings.cpp`
  - Exposes class `fftqspa.BCJRQSPA`.
  - Methods:
    - `info_bits_len()`
    - `code_bits_len()`
    - `encoder4bibo(info_bits)`
    - `decode4bibo(rr_bits_prob)`

- `python/setup.py`
  - Builds extension module `fftqspa` using MSVC flags.

- `python/main.py`
  - Example simulation script for DNA channel integration.
  - Uses:
    - `DNAChannel(...)` from `Inner_Code_DNA_Channel_Simulation.py`
    - `rr_bits_prob = 1 - voting_scores` then clipping.

- `python/Inner_Code_DNA_Channel_Simulation.py`
  - End-to-end inner-code + DNA-storage-channel pipeline.
  - The DNA channel model includes:
    - Synthesis stage with dropout and bias.
    - PCR amplification with bias and stochasticity.
    - Decay stage with dropout.
    - Sampling stage with random sampling of reads.
    - Sequencing stage with substitution errors.
  - Some of the Channel parameters are currently hardcoded in `main.py` and passed to `DNAChannel(...)`. Other Channel are listed in `python/Model/config.py`.
  - Calls MATLAB BCH encoder/decoder scripts under `python/Encode/`.

---

## Key Runtime Interfaces

### Python codec API contract

`codec = fftqspa.BCJRQSPA(parity_filename, max_iteration, mapping_filename)`

- `info_bits_len()`
  - Returns expected length of info bit vector (1-D).
  - Internally: `get_block_length() * m_degree`.

- `code_bits_len()`
  - Returns expected length of code bit vector (1-D).
  - Internally: `get_total_length() * m_degree`.

- `encoder4bibo(info_bits)`
  - Input: 1-D int array, values should be 0/1, length `info_bits_len()`.
  - Output: 1-D int array, length `code_bits_len()`, values 0/1.

- `decode4bibo(rr_bits_prob)`
  - Input: 1-D float array, length `code_bits_len()`.
  - Each entry means **P(bit=0)**.
  - Output: `(decoded_bits, iter)` where decoded bits are 0/1, length `code_bits_len()`.

### Decoder probability convention (critical)
- The decoder expects **probability domain**, not LLR.
- Input must be `P(bit=0)` per binary bit location.
- Avoid exact 0/1 probabilities; clip with epsilon (e.g., `1e-5`) to prevent numerical issues.

---

## LDPC Data File Expectations

The parity file loaded by `BCJRQSPA::Malloc(...)` is expected to include:
- `row_number`
- `col_number`
- `rank`
- `ary`
- sparse row-wise entries for parity matrix
- sparse row-wise entries for encoding/systematic matrix

`python/generate_binary_ldpc_dat.py` can convert a degree-3 text format into this `.dat` structure for binary (`q_ary=2`) cases.

---

## Build/Environment Notes

### C++ / Windows
- Existing VS Code task compiles active C++ file using `cl.exe`.
- Python extension build is configured in `python/setup.py` and compiles multiple C++ source files.

### Python dependencies
- See `python/requirements.txt`.
- Core required for extension usage and script flow:
  - `numpy`
  - `pybind11`
  - plotting/science stack as needed by scripts.

### MATLAB dependency
- DNA channel pipeline currently depends on MATLAB Engine (`matlab.engine`) and `.m` BCH scripts.
- If MATLAB is unavailable, DNA pipeline in `Inner_Code_DNA_Channel_Simulation.py` will not run as-is.

---

## Coding/Modification Rules for LLMs

1. Preserve old C++ style and memory model.
   - Manual `new[]/delete[]` is widely used.
   - Do not refactor broadly unless explicitly requested.

2. Keep algorithm semantics unchanged unless requested.
   - Especially FFT-QSPA message update order and normalization logic.

3. Validate shape/length contracts at Python boundary.
   - Binding currently enforces 1-D arrays and exact lengths; maintain this behavior.

4. Respect bit-order assumptions.
   - C++ packs/unpacks q-ary symbols using low-bit-first shifts (`(symbol >> ii) & 1`).

5. Maintain probability convention.
   - Any channel model should output or convert to `P(bit=0)` for `decode4bibo`.

6. Keep outputs compatible with existing result processing.
   - `results/*.csv` formats should remain stable unless requested.

---

## Common Pitfalls

- Mismatch between `q_ary` and mapping file signal count.
- Passing LLR into decoder instead of `P(bit=0)`.
- Providing 2-D arrays to bindings that require 1-D per codeword.
- Forgetting clipping of probabilities near 0/1.
- Assuming DNA channel scripts run without MATLAB Engine setup.
- Confusion between codeword axis and frame axis in numpy arrays (`(n_code, k_2)` pattern in scripts).

---

## Quick Start (Python simulation path)

1. Build/install extension from `python/`:
   - `python setup.py build_ext --inplace`
2. Run `python/main.py`.
3. Ensure parity/mapping files referenced in script exist in expected relative paths.
4. Ensure MATLAB + engine is available if DNA channel path is enabled.

---

## Important Project Assumptions (Current)

- Information bits are in the **tail** of systematic codeword layout in `python/main.py` post-processing (`sys_start = n_code - n_info`).
- `k_2` in `python/main.py` is used as simulation batch/frame count in loops, while also named like a data length in DNA channel code; treat with care and verify intent before changing.
- `innerRedundancy` may be partly overridden by fixed split (`r_1, r_2 = 15, 99`) in current DNA channel script.

---

## If You Are an LLM Working on This Repo

When asked to change this project, prioritize:
1. Keeping parity/mapping and array-shape contracts correct.
2. Preserving decoder input semantics (`P(bit=0)`).
3. Making minimal, targeted edits.
4. Calling out ambiguity explicitly (especially around `k_2`, redundancy allocation, and DNA-channel parameter meaning) before changing algorithm behavior.
