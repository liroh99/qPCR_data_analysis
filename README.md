# qPCR Relative Expression Analysis

## Overview

This project automates the analysis of qPCR (quantitative Polymerase Chain Reaction) data, which is often used in research to quantify gene expression levels. By automating the data extraction and analysis process, this project aims to save time and reduce potential errors associated with manual data processing. The program currently handles relative expression analysis (comparing treatment and control) but does not create standard curves.

## Scientific Background

Quantitative PCR (qPCR) is a laboratory technique used to measure the quantity of a specific RNA or DNA sequence in a sample. It is commonly used to quantify gene expression levels by measuring the cycle threshold (CT) values during PCR amplification.

This is how a standard qPCR curve looks like:
![alt text](https://www.bio-rad.com/sites/default/files/2022-03/identifying-cq-value-qpcr-high-res.jpg)

This qPCR graph illustrates the amplification of a DNA sample over a series of PCR cycles. The X-axis represents the number of cycles, while the Y-axis shows fluorescence intensity, which correlates with the amount of DNA produced. The baseline phase indicates the initial cycles with low fluorescence. The threshold line is set above the baseline to distinguish significant amplification. The Cq value (quantification cycle) is the cycle number where the sample's fluorescence (blue curve) surpasses the threshold, reflecting the point of detectable DNA amplification. The "No template" control shows no amplification, indicating no contamination in the reaction.


For more information about the scientific background, you can watch this YouTube video: [Analyzing Quantitative PCR Data](https://www.youtube.com/watch?v=y8tHiH0BzGY)

### Key Concepts

- **CT Value (Threshold Cycle)**: The cycle number at which the fluorescence signal of the PCR product crosses the detection threshold. Lower CT values indicate higher amounts of target nucleic acid.
- **Reference Gene**: A housekeeping gene used as a control for normalization of gene expression data.
- **Treated Groups**: Samples that have undergone a specific treatment or condition.
- **Water Template**: A negative control sample with water instead of a template to check for contamination.
- **Differential Expression Analysis**: A method to compare gene expression levels between different conditions (e.g., treated vs. control).



## Project Workflow

### Step 1: User Input

The program offers two ways to provide input:

1. **Graphical User Interface (GUI)**: When run without a config file, a GUI prompts for the following information:
   - Reference Gene
   - Control Sample Names (NTC, comma-separated)
   - Water Template Sample Name
   - Number of Treated Groups
   - Replicates (2 for duplicates, 3 for triplicates)
   - For each treated group:
     - Name of the treated group
     - Sample names for that group (comma-separated)

2. **Configuration File**: Alternatively, you can provide a YAML configuration file with the above information.

### Step 2: Data Processing

The program performs the following steps:

1. Loads data from the qPCR Excel file
2. Extracts relevant data (sample name, target name, CT values)
3. Calculates mean CT for each sample and target combination
4. Performs sanity checks (e.g., water template contamination)
5. Calculates ∆Cq, ∆Cq Expression, ∆∆Cq Expression, and % KD
6. Generates output table and plot

### Step 3: Output

The program generates two main outputs:

1. **CSV Table**: A CSV file containing the following columns:
   - Treatment
   - Sample
   - Target Gene
   - Cq Reference Gene
   - Cq Target
   - ∆Cq
   - ∆Cq Expression
   - Mean ∆Cq Expression
   - ∆Cq Expression stdev
   - ∆∆Cq Expression
   - ∆∆Cq Expression stdev
   - % KD

2. **Plot**: A bar plot showing the relative gene expression across treatments. The plot is saved as a JPEG file and can optionally be displayed on screen.

## Usage

To run the program, use the following command in the terminal:

```
python qPCR_data_analysis_with_errorbars --input_path <input_excel_file_path> --output_path <output_file_path> [options]
```

### Command-line Arguments

- `--input_path`, `-i`: Path to the qPCR raw data Excel file (required)
- `--output_path`, `-o`: Path for the output files (required)
- `--skiplines`: Number of lines to skip in the Excel file to reach column names (default: 47)
- `--config_path`, `-config`: Path to a YAML configuration file (optional, overrides GUI)
- `--show_plot`, `-show`: Flag to display the plot on screen (optional)

## Dependencies

To install the necessary dependencies, run:

```
pip install -r requirements.txt
```

Required packages include:
- pandas
- tkinter
- matplotlib
- numpy
- PyYAML

## Running the Tests

To run the tests, use:

```
pytest
```

## Note

This program is designed for relative expression analysis and does not support standard curve creation. Always verify the results and consult with a qPCR expert if you're unsure about the interpretation of the data.

This project was originally implemented as part of the Python programming course at the Weizmann Institute of Science taught by Gabor Szabo


