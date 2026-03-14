---
description: COPDS simulation platform context and coding rules for accurate LLM code generation
applyTo: "**/*"
---

# COPDS Project Instruction (Updated)

# Environment setup and migration requirements

This section is for a fresh machine (no Python installed) to ensure the project runs correctly.

## Required software

- Windows 64-bit.
- Python 3.12.x 64-bit from python.org (enable "Add Python to PATH" during install).
- Microsoft C++ Build Tools (MSVC v143) for compiling the pybind11 extension.

## Setup steps (PowerShell)

1. Create and activate a virtual environment in the repo root:
	- `python -m venv .venv`
	- `.\.venv\Scripts\Activate.ps1`
2. Install Python dependencies:
	- `python -m pip install --upgrade pip`
	- `python -m pip install -r requirements.txt`
3. Build the FFTQSPA extension in place:
	- `Set-Location Encode\Nonbinary`
	- `python setup.py build_ext --inplace`
	- Verify the output is `Encode/Nonbinary/fftqspa.cp312-win_amd64.pyd` for Python 3.12.
4. Return to repo root and run the app:
	- `Set-Location ../..`
	- `streamlit run Webapp.py`

## Notes and checks

- If the extension import fails, rebuild under `Encode/Nonbinary/` using the active venv.
- `pybind11` and `setuptools` are required for building the extension (included in `requirements.txt`).
- Use `Analysis/e2e_run_default_jpg.py` to validate an end-to-end run after setup.

## 1) Project purpose

COPDS is a Streamlit demo platform for DNA data storage simulation:
- encode image files into DNA strands,
- simulate DNA storage/sequencing channel effects,
- decode and reconstruct the original image.

Priority order:
1. End-to-end correctness.
2. Compatibility with existing artifacts and config structure.
3. Stable Streamlit behavior.
4. Readability and maintainability.

This is a demo/research workflow, not a production service.

## 2) Scope

### In scope
- Image files only: `jpg`, `jpeg`, `png`.
- Streamlit encode/restore workflow.
- Outer 16-ary LDPC + inner BCH pipeline.
- DNA channel simulation: synthesis, decay, PCR, sequencing.
- Reproducible intermediate files in project folders.

### Out of scope unless explicitly requested
- Arbitrary binary file support.
- Large UI redesign.
- Major architecture rewrite.

## 3) Runtime and environment facts

- Main language: Python.
- UI: Streamlit.
- Core backend: `Encode/Nonbinary/fftqspa` pybind11/C++ extension.
- Verified working environment in this repo: Python 3.12 virtual environment.
- For Python 3.12, the compiled extension should be `Encode/Nonbinary/fftqspa.cp312-win_amd64.pyd`.
- If the extension import fails, rebuild under `Encode/Nonbinary/` using the active environment.
- Required Python build/runtime dependencies discovered in this session: `setuptools`, `pybind11`.

Important:
- Do not hard-code a `cp311` extension path in new logic.
- Prefer version-aware module path generation when metadata must be written to config.

## 4) Main files and roles

- `Webapp.py`: Streamlit entry point and encode/restore orchestration.
- `Helper_Functions.py`: segmentation, DNA/binary transforms, BCH/LDPC wrappers, voting helpers.
- `Encode/Nonbinary/`: FFTQSPA codec sources, parity files, mapping files, compiled extension.
- `Model/Model.py`: DNA channel simulation.
- `Model/config.py`: default channel parameters and substitution matrices.
- `config/config.json`: persistent metadata for stored files and current ECC settings.
- `DNA_Library/`: `{file_id}_in.dna`, `{file_id}_out.dna`.
- `Mid_data/`: intermediate BCH/LDPC text files.
- `Restored_files/`: restored output files.
- `Analysis/e2e_run_default_jpg.py`: end-to-end validation script added during this session.

## 5) Canonical pipeline

### Encode path
1. User selects mode `Encode & Store in DNA`.
2. User uploads images and submits.
3. Each image is split by raw file bytes into chunks of `chunk_size / 8` bytes.
4. Chunks are converted to bit strings.
5. Chunks are grouped into pools of size `k_0`; the final pool is padded with random bit chunks.
6. Pools are transposed for inter-oligo LDPC encoding.
7. Outer LDPC encodes each row.
8. BCH encodes both index part and LDPC payload part.
9. The full BCH codeword (`index + data`) is converted to DNA.
10. Pre-channel DNA is saved to `DNA_Library/{file_id}_in.dna`.
11. DNA channel simulation runs.
12. Simulated DNA is saved to `DNA_Library/{file_id}_out.dna`.
13. Metadata is merged into `config/config.json`.

### Restore path
1. User selects mode `Restore from DNA`.
2. User selects stored file(s) and submits.
3. App loads `config/config.json` and `{file_id}_out.dna`.
4. DNA is converted back into binary arrays: `[ids, info, counts]`.
5. BCH decode + voting (`bch_decode_and_vote_cpp`) produces LDPC input probabilities.
6. LDPC decode reconstructs original chunk bits.
7. Chunks are converted back to bytes.
8. Padding chunks and final-byte padding are removed.
9. Restored file is written to `Restored_files/`.

## 6) Critical data contracts

### `config/config.json`
Must preserve this structure:
- `file_num`
- `files`
	- `{file_id}`
		- `file_name`
		- `file_type`
		- `suffix`
		- `file_size`
		- `segmentation.chunk_num`
		- `segmentation.block_num`
		- `segmentation.pad_size`
- `ECC.outer`
	- `n`, `k`, `h_matrix_path`, `ldpc_encoder_exe`, `parity_path`, `mapping_path`, `max_iter`, `code_ary`
- `ECC.inner1`
	- `n`, `k`
- `ECC.inner2`
	- `n`, `k`
- `ECC.BCH_codec_exe`
- `DNA2Binary.is_padding`
- `DNA2Binary.padding`

### File naming
- Input DNA pool: `{file_id}_in.dna`
- Simulated DNA pool: `{file_id}_out.dna`
- Mid data files use `{file_id}` with LDPC/BCH stage suffixes.

### DNA pool format
- Current valid encode output uses full BCH codewords, not payload-only DNA.
- `*_out.dna` is JSON where each key is `chunk_{block_id}` and each value is a dictionary:
	- key: DNA sequence string
	- value: count

### `_in.dna` writing rule
- Write a single valid JSON object, not multiple JSON fragments appended to one file.

## 7) Fixed parameters in current UI flow

These are effectively fixed in the Streamlit app:
- Outer LDPC: $(n_0, k_0) = (2048, 1024)$
- Inner BCH index: $(n_1, k_1) = (29, 14)$
- Inner BCH payload: $(n_2, k_2) = (419, 320)$ for current standard pipeline
- Parity file: `Encode/Nonbinary/Parity_files_2048/512_256_16aryCode17.dat`
- Mapping file: `Encode/Nonbinary/Mapping_files/SignalSet_BPSK-4.txt`

Decoder expectation:
- LDPC decode consumes $P(bit=0)$ probabilities.
- Probabilities should be clipped away from exactly 0 and 1.

## 8) Verified fixes from this session

These points reflect the current known-good behavior and should not be regressed:

1. **Streamlit startup works** under the current `.venv`.
2. **`file_uploader` state bug fixed**: do not mutate `st.session_state["file_uploader"]` after widget creation.
3. **Encode/restore actions must be submit-gated**: do not run full encode/decode on every rerun.
4. **Config merge behavior fixed**: do not overwrite `config/config.json` before merging file records.
5. **DNA conversion root-cause fixed**: DNA must be generated from full BCH output (`index + data`), not payload-only BCH text.
6. **LDPC decode output orientation fixed**: do not transpose decoded rows before writing restore-stage output.
7. **DNA count aggregation fixed**: repeated simulated DNA sequences must accumulate counts.
8. **Restore multi-file chunk contamination fixed**: reset decoded chunk buffers per file.
9. **Missing file/missing chunk guards added**.

## 9) Legacy compatibility warning

Some historical DNA pools in `DNA_Library/` are not decodable with the corrected pipeline because they were generated by the old broken encode path.

Known cause:
- old DNA pools may contain only BCH payload bits and miss BCH index bits.

Rule:
- if a historical pool is detected as legacy payload-only format, do not silently decode it;
- show a clear error and require re-encoding with the fixed pipeline.

Concrete example from this session:
- historical `Fruits.jpg` pool is legacy/incompatible and must be re-encoded.

## 10) End-to-end validation result from this session

Verified successfully:
- `Test_Files/default.jpg` was re-encoded and restored end-to-end.
- Restored file matched the source exactly by SHA256.

Validation script:
- `Analysis/e2e_run_default_jpg.py`

Use that script as the fastest known-good non-UI regression check.

## 11) Coding rules for future edits

1. Prefer small, surgical changes.
2. Preserve config keys, file names, and intermediate file contracts.
3. Use project-relative paths.
4. Do not introduce hard-coded external absolute paths.
5. Keep Streamlit flow simple and stable.
6. When changing encode/decode logic, validate with an end-to-end run.
7. Keep FFTQSPA integration working; do not break `Encode/Nonbinary/` import/build.
8. When touching restoration logic, consider both current correct pools and legacy broken pools.

## 12) Minimal validation checklist

After encode/decode changes, verify:
- Streamlit app launches.
- FFTQSPA imports correctly.
- Encode produces valid `*_in.dna` and `*_out.dna`.
- Restore reconstructs at least one newly encoded image successfully.
- `config/config.json` remains parseable and contract-compatible.
- `Analysis/e2e_run_default_jpg.py` still succeeds.

If intent is unclear, default to minimal compatible changes.