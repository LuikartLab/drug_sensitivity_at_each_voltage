# Provenance / AI-assisted development

## Summary
This repository contains analysis scripts for drug sensitivity calculations from I/V curve CSV files.

The script `DrugSensitivityAtEachVoltage.py` was developed through an iterative, tool-assisted workflow: I (the repository author) used ChatGPT as a coding assistant to draft and refine the script, then I tested it locally, reported results/errors, and requested modifications until outputs matched the intended analysis and visualization requirements.

ChatGPT did not run the electrophysiology experiments or validate biological interpretations. I am responsible for verifying correctness, maintaining the code, and interpreting results.

## Script covered
- `DrugSensitivityAtEachVoltage.py`

## Development process (high level)
Key iterations included:

1. **Initial AUC-based analysis (baseline)**  
   - Started from an earlier script that computed AUC across a voltage range for WT/dKO control vs drug.

2. **Switch to per-voltage current-based analysis**  
   - Updated the computation to use current at each voltage step rather than AUC.
   - Output expanded to “one row per voltage per file” including mean/SEM, block%, and sensitivity ratio.

3. **Heatmap generation for ImageJ**  
   - Added 8-bit heatmap export where x-axis is voltage and y-axis is drug/file.
   - Added scaling to 8-bit and increased cell size (e.g., 100×100 pixels per cell) for clearer visualization.

4. **Labeled visualization outputs**  
   - Added a labeled PNG with voltage/drug axes and a colorbar.

5. **Custom fire/ice LUT**  
   - Implemented a diverging LUT with black near zero, icy colors for negative values, and hot colors for positive values.
   - Adjusted “fire” to use darker reds through bright orange/yellow (no white endpoint).
   - Wrote an ImageJ-compatible `.lut` file.

6. **Bug fix: sensitivity ratio sign artifact**  
   - Discovered that small negative block values (apparent facilitation/noise) could yield negative WT/dKO ratios.
   - Updated sensitivity ratio to use magnitude: `abs(WT block) / abs(dKO block)` to reflect “sensitivity” as effect size.

## Evidence / record
- Development was performed via an interactive ChatGPT session.  
- If needed, a transcript excerpt or export can be added under `docs/chatgpt-transcript.md` (optional).

## How to cite / acknowledge
Recommended language:
> “This code was developed by the repository author with iterative assistance from ChatGPT (OpenAI) for drafting and refactoring.”
