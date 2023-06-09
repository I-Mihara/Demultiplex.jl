# Demultiplex
Demultiplex is a Julia package that provides functions for demultiplexing reads based on barcodes. The package includes two functions: "demultiplex" and "BCSemiglobalAlignmentScore".
## Installation
To install the Demultiplex package, open the Julia REPL and run the following command:
```Julia
using Pkg;Pkg.add("Demultiplex")
```
## Functions
### demultiplex
The "demultiplex" function takes four input arguments: "input"(R1 sequence), "input2"(R2 sequence), "bc_tsv", and "output_dir", as well as an optional fifth argument "maximum_errorrate". The function compares the reverse complementary stand converted barcode sequences in .tsv file  to a R1 sequences in input fastq file. It outputs files containing the R2 reads that correspond the reference barcodes.
```julia
function demultiplex(input::String,input2::String,bc_tsv::String,output_dir::String,maximum_errorrate=0.22::Float64)
```

#### Arguments
* `input::String`: the path to the input R1 fastq file
* `input2::String`: the path to the input R2 fastq file
* `bc_tsv::String`: the path to the reference barcode file
* `output_dir::String`: the path to the output directory where the demultiplexed files will be written.
* `maximum_errorate=0.22::Float64`: an optional argument that sets the maximum error rate for barcode matching. The default value is 0.22

### BCSemiglobalAlignmentScore
The "bcSemiglobalAlignmentScore" function takes threee input arguments: "ref", "query", and "maximum_errorrate". The function calculates a semi-global alignment score between the two input sequences,"ref" and "query". It then computes the error rate of the alignment and adjusts the score based on the number of matches, mismatches, and insertions. If the error rate is greater than the maximum allowed error rate, the score is set to zero. The function returns the final alignment score.

```julia
function BCSemiglobalAlignmentScore(ref::String,query::String,maximum_errorrate::Float64)
```
#### Arguments
* `ref`: the reference sequence.
* `query`: the query sequence to be aligned to the reference.
* `maximum_errorrate`: the maximum allowed rate for barcode matching.

## Usage
Here is an example of how to use the `demultiplex()` function
```julia
using Demultiplex

input = "path/to/input.fastq"
input2 = "path/to/input2.fastq"
bc_tsv = "path/to/barcodes.tsv"
output_dir = "path/to/output/directory"
maximum_errorrate = 0.22

demultiplex(input, input2, bc_tsv, output_dir, maximum_errorrate)
```
To use the function in the terminal, install the "Demultiplex" package by running the following command first
```
julia -e 'using Pkg; Pkg.add("Demultiplex")'
```
Once installed, run the `demultiplex()` function with your desired inputs using this command:
```
julia -e 'using Demultiplex; demultiplex(input, input2, bc_tsv, output_dir, maximum_errorrate)'
```

## License
This package is licensed under the MIT License.

