# qPCR Project

## Overview

This project automates the analysis of qPCR data, which is often used in research to quantify gene expression levels. By automating the data extraction and analysis process, this project aims to save time and reduce potential errors associated with manual data processing.

## Scientific Background

Quantitative PCR (qPCR) is a laboratory technique used to measure the quantity of a specific RNA or DNA sequence in a sample. It is commonly used to quantify gene expression levels by measuring the cycle threshold (CT) values during PCR amplification.

### Key Concepts

- **CT Value (Threshold Cycle)**: The cycle number at which the fluorescence signal of the PCR product crosses the detection threshold. Lower CT values indicate higher amounts of target nucleic acid.
- **Reference Gene**: A housekeeping gene used as a control for normalization of gene expression data.
- **Treated Groups**: Samples that have undergone a specific treatment or condition.
- **Water Template**: A negative control sample with water instead of a template to check for contamination.
- **Differential Expression Analysis**: A method to compare gene expression levels between different conditions (e.g., treated vs. control).

## Project Workflow

### Step 1: User Input via GUI

When you run the program, a graphical user interface (GUI) will prompt you to provide the following information:

1. **Reference Gene**: The name of the housekeeping gene used for normalization.
2. **Control Sample Name (NTC)**: The name of the non-treated control sample.
3. **Water Template Sample Name**: The name of the water template sample.
4. **Number of Treated Groups**: The number of treated groups.
5. **Treated Groups Information**: For each treated group, you will be asked to provide:
   - The name of the treated group.
   - The sample names corresponding to that treated group (comma-separated).

### Step 2: Running the Program

To run the program, use the following command in the terminal:
```sh
python qpcr_analysis.py --input <input_excel_file_path> --output <output_excel_file_name>


### Step 3: Running the Program

The program performs the following steps:

1. **Load Data**: Reads the qPCR Excel file and focuses on the "Results" sheet.
2. **Extract Relevant Data**: Extracts the sample name, target name, CT values, and CT SD values from the "Results" sheet and creates a DataFrame.
3. **Calculate Mean and Standard Deviation**: Calculates the mean and standard deviation of CT values for duplicates/triplicates.
4. **Check Water Template**: Verifies if the water template is clean (CT > 35). If not, it displays a warning.
5. **Check CT SD**: Ensures that the CT SD is within the acceptable range (CT SD ≤ 0.5). If any sample exceeds this threshold, it displays a warning.
6. **Calculate ∆Cq**: Calculates the difference between the CT values of the target gene and the reference gene for each sample.
7. **Calculate Expression (2^(-∆Cq))**: Computes the expression value for each sample.
8. **Calculate Mean and Standard Deviation of Expression**: Calculates the mean and standard deviation of the expression values for duplicates/triplicates.
9. **Calculate ∆∆Cq**: Computes the ∆∆Cq expression for treated vs. control.
10. **Calculate ∆∆Cq Standard Deviation**: Calculates the standard deviation for ∆∆Cq expression.
11. **Calculate % KD**: Computes the knockdown percentage.
12. **Generate Output**: Creates a table similar to the provided example and generates a graph.

### Output

- **Table**: An Excel file containing a table with the following columns:
  - Treatment
  - Cq reference gene
  - Cq target
  - ∆Cq
  - ∆Cq expression
  - Mean ∆Cq expression
  - ∆Cq expression stdev
  - ∆∆Cq expression
  - ∆∆Cq expression stdev
  - % KD

- **Graph**: A graph where the NTC is set to 1, and the treated groups are represented by the ∆∆Cq expression values. This graph shows the differential expression of the gene of interest due to the treatment.

## Dependencies

To install the necessary dependencies, run:
```sh
pip install -r requirements.txt
