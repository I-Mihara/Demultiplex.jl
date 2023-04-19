# DemultiplexJulia
DemultiplexJulia is a Julia package that provides functions for demultiplexing reads based on barcodes. The package includes two functions: "Demultiplexer" and "BBCSemiglobalAlignmentScore".
## Installation
To install the DemultiplexJulia package, open the Julia REPL and run the following command:
```Julia
using Pkg;Pkg.add DemultiplexJulia
```
## Functions
### Demultiplexer
The "Demultiplexer" function takes four input arguments:"input","input2","bbc_tsv",and"output_dir",and an optional fifth argument "maximum_errorrate". The function reads in fastq files compares the barcode sequences to a reference file, and outputs files containing the reads that match the reference barcodes.
```julia
function Demultiplexer(input::String,input2::String,bbc_tsv::String,output_dir::String,maximum_errorrate=0.22::Float64)
```

#### Arguments
* `input::String`:the path to the input R1 fastq file
* `input2::String`:the path to the input R2 fastq file
* `bbc_tsv::String`:the path to the reference barcode file
* `output_dir::String`: the path to the output directory where the demultiplexed files will be written.
* `maximum_errorate=0.22::Float64`:an optional argument that sets the maximum error rate for barcode matching. The default value is 0.22

### BBCSemiglobalAlignmentScore
The "BBCSemiglobalAlignmentScore"function takes threee input arguments:"ref","query",and"maximum_errorrate".The function calculates a semi-global alignment score between the two input sequences,"ref" and "query", using an affine gap score model. It then computes the error rate of the alignment and adjusts the score based on the number of matches, mismatches, and insertions. If the error rate is greater than the maximum allowed error rate, the score is set to zero. The function returns the final alignment score.

```julia
function BBCSemiglobalAlignmentScore(ref::String,query::String,maximum_errorrate::Float64)
```
#### arguments
* `ref`:the reference sequence.
* `query`:the query sequence to be aligned to the reference.
* `maximum_errorrate`:the maximum allowed rate for barcode matching.

## Usage
here is an example of how to use the "Demultiplexer" function
```julia
using DemultiplexJulia

input = "path/to/input.fastq"
input2 = "path/to/input2.fastq"
bbc_tsv = "path/to/barcodes.tsv"
output_dir = "path/to/output/directory"
maximum_errorrate = 0.2

Demultiplexer(input, input2, bbc_tsv, output_dir, maximum_errorrate)
```
## License
???

