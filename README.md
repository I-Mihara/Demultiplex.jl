# Demux.jl

## Overview
Demux.jl is a Julia package designed for demultiplexing reads based on barcodes. Given a set of sequencing reads in FASTQ format and a reference barcode in TSV format, each read is assigned to a FASTQ file corresponding to its barcode. This barcode assignment process is designed to be robust to barcode mutations and allows you to adjust the permitted level of mutation through parameters.

### Program features
* Fast and accurate semi-global alignment
* Robust to barcode mutations
* No restrictions on barcode size or position
* Usable as julia package
* Supports parallel computing

### References

## Installation

### Step 0: Install Julia
If you haven't installed Julia yet, you can download it from the [official Julia website](https://julialang.org/downloads/). Follow the instructions there for your operating system.

### Option 1: Using Julia's Pkg REPL Mode
1. Open the Julia REPL (by typing `julia` in the terminal).
2. Press `]` to enter the Pkg REPL mode.
3. Run the following command to install `Demux.jl`:
    ```julia
    add Demux
    ```
4. Press the `Backspace` key to exit the Pkg mode.

### Option 2: Using the Command Line
1. Open a terminal or command prompt.
2. Run Julia and execute the installation command directly:
    ```bash
    julia -e 'using Pkg; Pkg.add("Demux")'
    ```

Both methods will install `Demux.jl` and its dependencies, making it ready for use in your projects.

### Dependencies
* Julia >= 1.10.5 (includes the `Distributed` standard library)
* DataFrames >= 1.7.0
* CSV >= 0.10.14

## Basic Usage

The primary function of this package is `execute_demultiplexing()`. It can classify sequences in FASTQ file using barcodes from reference file. Usage is as follows:
```Julia
execute_demultiplexing(FASTQ_file, barcode_file, output_directory)
```

### Input

#### FASTQ File
* There is no restriction on the sequence length in the FASTQ file.
* The function can take one or two FASTQ files as input. In the case of using two FASTQ files (R1 and R2), the command can be executed as follows:
```julia
execute_demultiplexing(FASTQ_R1, FASTQ_R2, barcode_file, output_directory)
```
When using two FASTQ files, sequences in the R2 file are classified by calculating similarity scores from the R1 sequences and barcodes in the reference file.

#### Barcode Reference File
* The reference file is expected to be a TSV file containing the following columns: `ID`, `Full_seq`, `Full_annotation`, as shown below:
```
ID  Full_seq	Full_annotation
001-barcode ACAGACUACAAA LLLBBBBBBB33
```
* In the `Full_seq` column, the region specified as `B` in the `Full_annotation` column is considerd as the barcode.

### Output

* All output files will be saved in the specified `output_directory`.
* The names of the output files are based on the `ID` values in the barcode reference file. For example, if the reference file contains IDs such as `001` and `002`, the resulting output files will be named `001.fastq`, `002.fastq`, and so on.
* Sequences that do not match any barcode in the reference file are saved in `unknown.fastq`. Sequences that have ambiguous classification (i.e., they match multiple barcodes with similar scores) are saved in `ambiguous_classification.fastq`.

## Options

The `execute_demultiplexing` function provides several optional parameters to control the demultiplexing process:

- **`max_error_rate::Float64`** (default: `0.22`): 
  - The maximum allowable error rate when matching sequences to barcodes. A lower value increases stringency, allowing fewer mismatches and indels. Sequences that exceed this error rate threshold and do not match any barcode will be saved in `unknown.fastq`.
  
- **`min_delta::Float64`** (default: `0`): 
  - The minimum difference in similarity scores required to classify a sequence. This helps in resolving cases where a sequence may match multiple barcodes. A larger value increases the threshold, reducing ambiguous classifications. Sequences that match multiple barcodes and do not meet this minimum difference will be saved in `ambiguous_classification.fastq`.
  
- **`mismatch::Int`** (default: `1`): 
  - The penalty score for mismatches during sequence alignment. Increasing this value makes the alignment process more stringent by penalizing mismatches more heavily.
  
- **`indel::Int`** (default: `1`): 
  - The penalty score for insertions and deletions (indels) during sequence alignment. A higher value imposes a stricter penalty for gaps in the alignment.

- **`classify_both::Bool`** (default: `false`): 
  - If set to `true`, the function will classify both R1 and R2 sequences and output separate files for each. Otherwise, it classifies only R2 sequences by default.

- **`bc_rev::Bool`** (default: `false`): 
  - If set to `true`, the barcodes in the reference file are reversed before comparison. This is useful when the barcode sequences are provided in reverse orientation.

### Behavior
- **Output Directory**: 
  - The function will throw an error if the specified `output_dir` already exists to prevent overwriting. A new directory is created to store the output files.

- **Parallel Processing**: 
  - The function automatically detects the number of available workers (`nworkers()`) for parallel processing. When there is more than one worker, the function divides the input FASTQ files into parts for concurrent processing, accelerating the demultiplexing process.

### Usage Example
```julia
execute_demultiplexing(FASTQ_R1, FASTQ_R2, barcode_file, output_directory, max_error_rate=0.2, min_delta=0.1, mismatch=1, indel=2, classify_both=true, bc_rev=true)
```
