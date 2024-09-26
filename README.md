# Demux.jl

## Overview
Demux.jl is a Julia package designed for demultiplexing reads based on barcodes. Given a set of sequencing reads in FASTQ format and a reference barcode in TSV format, each read is assigned to a FASTQ file corresponding to its barcode. This barcode assignment process is designed to be robust to barcode mutations and allows you to adjust the permitted level of mutation through parameters.

### Program features
* Fast and accurate semi-global alignment
* Robust to barcode mutations
* No restriction s on barcode size or position
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

