---
description: COPDS simulation platform context and coding rules for accurate LLM code generation
applyTo: "**/*"
---

# COPDS Project Instruction (Simulation Platform with UI)

## 1) Project mission

This project is an interactive **demonstration platform** for DNA data storage simulation.

Core goal:
- Encode uploaded images into DNA strands.
- Simulate DNA storage channel errors.
- Decode and restore images from simulated DNA data.

This is not a production backend service; prioritize:
- Correct end-to-end behavior.
- Stable Streamlit UX.
- Reproducible intermediate artifacts.

## 2) Current scope and product position

### In scope
- Image files only (`jpg`, `jpeg`, `png`).
- Streamlit-based encode/restore workflow.
- Two-layer BCH + outer LDPC pipeline.
- DNA channel simulation (synthesis, decay, PCR, sequencing).

### Out of scope (unless user explicitly asks)
- General arbitrary binary file support.
- Major UI redesign.
- Full architecture rewrite.

## 3) Technology and runtime constraints

- Language: Python (main orchestration/UI).
- UI framework: Streamlit.
- Numeric stack: NumPy, SciPy ecosystem.
- Visualization: Plotly / Matplotlib / Seaborn.
- Current codec backend: external executables under `Encode/` (e.g., LDPC/BCH).

Important roadmap direction:
- Keep current executable-based pipeline working.
- Prefer modular changes that make it easy to **gradually replace codec components with Python + C++ implementations**.
- Do not break existing interfaces while migrating.

## 4) Entry points and important modules

- `Webapp.py`: main Streamlit app, mode switch (`Encode & Store in DNA` / `Restore from DNA`).
- `Helper_Functions.py`: data transforms, segmentation, DNA conversion, voting, codec wrappers.
- `Model/Model.py`: DNA channel simulation modules.
- `Model/config.py`: default channel parameters and substitution transition matrices.
- `config/config.json`: persistent coding metadata for stored files.
- `DNA_Library/`: generated DNA pools (`*_in.dna`, `*_out.dna`).
- `Mid_data/`: intermediate files for BCH/LDPC encode/decode IO.
- `Restored_files/`: reconstructed output files.

## 5) End-to-end data flow (must stay consistent)

### Encode path
1. Upload one or more images from Streamlit.
2. Split image bytes into fixed-size chunks (`chunk_size / 8` bytes each).
3. Convert bytes to bit strings.
4. Group chunks into blocks of size `k_0`.
5. Outer LDPC encode per block.
6. Inner BCH encode (index + payload).
7. Convert binary codewords to DNA sequences.
8. Save pre-channel DNA pools to `DNA_Library/*_in.dna`.
9. Run channel simulation and save post-channel pools to `DNA_Library/*_out.dna`.
10. Persist metadata to `config/config.json`.

### Restore path
1. Load selected file metadata from `config/config.json`.
2. Load corresponding `*_out.dna` pools.
3. Convert DNA back to binary (split index/info/counts).
4. BCH decode index and payload.
5. Vote by index to build probability matrix.
6. LDPC decode.
7. Reassemble chunks and remove padding.
8. Write restored files to `Restored_files/` and display images.

## 6) Data contracts (critical)

### `config/config.json` minimum structure
Must preserve keys and nesting:
- `file_num`
- `files` (map: `file_id -> file metadata`)
  - `file_name`, `file_type`, `suffix`, `file_size`
  - `segmentation.chunk_num`, `segmentation.block_num`, `segmentation.pad_size`
- `ECC.outer` (`n`, `k`, `h_matrix_path`, `ldpc_encoder_exe`)
- `ECC.inner1` (`n`, `k`)
- `ECC.inner2` (`n`, `k`)
- `ECC.BCH_codec_exe`
- `DNA2Binary.is_padding`, `DNA2Binary.padding`

### DNA library file naming
- Input pool: `{file_id}_in.dna`
- Simulated output pool: `{file_id}_out.dna`

### Mid data naming (do not change lightly)
- LDPC encode/decode and BCH encode/decode files use `{file_id}` + block/chunk indices.
- Any migration must keep backward compatibility or include explicit conversion logic.

## 7) LLM coding rules for this repository

1. **Preserve pipeline compatibility first**.
	- Do not rename core config keys, file naming conventions, or directory contracts.

2. **Prefer small, surgical edits**.
	- Avoid broad refactors unless explicitly requested.

3. **Use path-safe logic**.
	- Prefer project-relative paths with `os.path.join`.
	- Avoid introducing new hard-coded absolute paths.

4. **Keep Streamlit interactions stable**.
	- Retain mode switch behavior and basic form flow.
	- Do not add unrelated UI complexity.

5. **Guard edge cases for robustness**.
	- Empty upload list.
	- Missing config entries.
	- Missing DNA files.
	- Padding handling in final chunk.

6. **Design for codec migration (Python + C++)**.
	- Isolate codec calls behind clear interfaces.
	- Keep adapter-like boundaries so executable backends can be swapped incrementally.

7. **Maintain existing output artifacts**.
	- Keep generated files in `DNA_Library/`, `Mid_data/`, and `Restored_files/` unless task requires otherwise.

## 8) Preferred implementation style

- Keep functions cohesive and single-purpose.
- Use explicit variable names (avoid one-letter names unless mathematically standard).
- Add docstrings for new non-trivial functions.
- Log/print concise progress at critical stages (encode, channel, decode).
- Keep bilingual comments only when necessary; prefer concise English technical wording.

## 9) Validation checklist for code changes

When changing encode/decode logic, verify:
- Streamlit app starts successfully.
- Encode path can produce `*_in.dna` and `*_out.dna`.
- Restore path can reconstruct at least one previously encoded image.
- `config/config.json` remains parseable and contract-compatible.
- No breaking change to BCH/LDPC IO file format unless explicitly planned.

## 10) Priority order for decision making

When requirements conflict, follow this order:
1. End-to-end correctness of encode/restore pipeline.
2. Compatibility with existing data/config/artifacts.
3. Stability of Streamlit user flow.
4. Readability and maintainability.
5. Performance optimization.

## 11) What the LLM should ask before large changes

Before major modifications, clarify:
- Whether change targets demo stability or migration toward Python/C++ codecs.
- Whether backward compatibility with existing `DNA_Library` and `config/config.json` is required.
- Whether scope is UI-only, algorithm-only, or both.

If user intent is unclear, default to minimal compatible changes.