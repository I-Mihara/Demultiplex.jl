module Demultiplex

export
    BCSemiglobalAlignmentScore,
    demultiplex

using BioAlignments, CSV, DataFrames

"""
    BCSemiglobalAlignmentScore(ref::String,query::String,maximum_errorrate=0.22::Float64) -> score::Int64 
Calculate the semiglobal alignment score between a reference sequence and a query sequence
using the Bioalignments.jl package and return score based on the number of mach, mismach and insertion.
If the error rate between the aligned seqences is greater than the maximum allowed error rate,
the function returns a score of 0.
"""
function BCSemiglobalAlignmentScore(ref::String, query::String, maximum_errorrate=0.22::Float64)
    problem = SemiGlobalAlignment()
    scoremodel = AffineGapScoreModel(match=0, mismatch=-1, gap_open=0, gap_extend=-1)
    result = pairalign(problem, query, ref, scoremodel)
    alignment_result = alignment(result)
    errorrate = -BioAlignments.score(result) / length(query)
    score = count_matches(alignment_result) - count_mismatches(alignment_result) - 2 * count_insertions(alignment_result)
    if errorrate > maximum_errorrate
        score = 0
    end
    return score
end

"""
    function demultiplex(input::String,input2::String,bc_tsv::String,output_dir::String,maximum_errorrate=0.22::Float64)
Demultiplex reads based on barcode sequences, with the option to set a maximum error rate for matching.
"""
function demultiplex(input::String, input2::String, bc_tsv::String, output_dir::String, maximum_errorrate=0.22::Float64)
    if isdir(output_dir)
        rm(output_dir, recursive=true)
    end
    mkdir(output_dir)
    bc_df = CSV.read(bc_tsv, DataFrame, delim="\t")
    for i in 1:nrow(bc_df)
    prefix_region = 1:findfirst('B', bc_df.Full_annotation[i])-1
    suffix_region = findlast('B', bc_df.Full_annotation[i])+1:length(bc_df.Full_annotation[i])
    prefix = SubString(bc_df.Full_seq[i], prefix_region)
    suffix = SubString(bc_df.Full_seq[i], suffix_region)
    bc_df.Full_seq[i] = replace(bc_df.Full_seq[i], prefix => "")
    bc_df.Full_seq[i] = replace(bc_df.Full_seq[i], suffix => "")
    end
    bc_df.Full_seq = replace.(bc_df.Full_seq, "U" => "T")
    bc_df.Full_seq = reverse.(bc_df.Full_seq)
    bc_df.Full_seq = replace.(bc_df.Full_seq, "A" => "T", "T" => "A", "G" => "C", "C" => "G")
    open(input, "r") do file
        mode = ""
        header = ""
        filename = ""
        seq = ""
        line3 = ""
        quality_score = ""
        open(input2, "r") do file2
            for line in eachline(file)
                R2 = readline(file2)
                if line[1] == '@' && mode == ""
                    header = R2
                    mode = "header"
                elseif mode == "header"
                    seq = line
                    mode = "line3"
                    maximum_bc_score = 0
                    maximum_bc_score_number = 0
                    for (i, row) in enumerate(eachrow(bc_df))
                        bc_score = BCSemiglobalAlignmentScore(seq, row.Full_seq, maximum_errorrate)
                        if maximum_bc_score < bc_score
                            maximum_bc_score = bc_score
                            maximum_bc_score_number = i
                        end
                    end
                    seq = R2
                    if maximum_bc_score_number == 0
                        filename = output_dir * "/trimmed-unknown.fastq"
                    else
                        filename = output_dir * "/trimmed-" * string(bc_df.ID[maximum_bc_score_number]) * ".fastq"
                    end
                elseif mode == "line3"
                    line3 = R2
                    mode = "quality_score"
                elseif mode == "quality_score"
                    quality_score = R2
                    open(filename, "a") do outputfile
                        write(outputfile, header * "\n" * seq * "\n" * line3 * "\n" * quality_score * "\n")
                    end
                    mode = ""
                end
            end
        end
    end
end
end
