module DemultiplexJulia

export 
BBCSemiglobalAlignmentScore,
Demultiplexer


#include("/home/ec2-user/BioAlignment/src/BioAlignments.jl")
#include("")
using BioAlignments , CSV , DataFrames#使う関数だけusingする？

"""
    BBCSemiglobalAlignmentScore(ref::String,query::String,maximum_errorrate=0.22::Float64) -> score::Int64 
Calculate the semiglobal alignment score between a reference sequence and a query sequence
using the Bioalignments.jl package and return score based on the number of mach, mismach and insertion.
If the error rate between the aligned seqences is greater than the maximum allowed error rate,
the function returns a score of 0.
"""
function BBCSemiglobalAlignmentScore(ref::String,query::String,maximum_errorrate=0.22::Float64)
    problem = SemiGlobalAlignment()
    scoremodel = AffineGapScoreModel(match=0,mismatch=-1,gap_open=0,gap_extend=-1)
    result=Bioalignments.pairalign(problem, query, ref, scoremodel)
    alignment_result = alignment(result)
    errorrate = - BioAlignments.score(result) / length(query)
    score = count_matches(alignment_result)-count_mismatches(alignment_result)-2*count_insertions(alignment_result)
    if errorrate > maximum_errorrate
        score = 0
    end
    return score
end

"""
    function Demultiplexer(input::String,input2::String,bbc_tsv::String,output_dir::String,maximum_errorrate=0.22::Float64)
Demultiplex reads based on barcode sequences, with the option to set a maximum error rate for matching.
"""
function Demultiplexer(input::String,input2::String,bbc_tsv::String,output_dir::String,maximum_errorrate=0.22::Float64)
    if isdir(output_dir)
        rm(output_dir,recursive=true) 
    end
    mkdir(output_dir)
    bbc_df = CSV.read(bbc_tsv,DataFrame,delim="\t")
    prefix_region = 1:findfirst('B',bbc_df.Full_annotation[1])-1
    suffix_region = findlast('B',bbc_df.Full_annotation[1])+1:length(bbc_df.Full_annotation[1])
    prefix = SubString(bbc_df.Full_seq[1],prefix_region)
    suffix = SubString(bbc_df.Full_seq[1],suffix_region)
    bbc_df.Full_seq = replace.(bbc_df.Full_seq,prefix=>"")
    bbc_df.Full_seq = replace.(bbc_df.Full_seq,suffix=>"")
    bbc_df.Full_seq = replace.(bbc_df.Full_seq,"U"=>"T")
    bbc_df.Full_seq = reverse.(bbc_df.Full_seq)
    bbc_df.Full_seq = replace.(bbc_df.Full_seq,"A"=>"T","T"=>"A","G"=>"C","C"=>"G")
    open(input,"r") do file
        mode = ""
        header = ""
        filename = ""
        seq = ""
        line3 = ""
        quality_score = ""
        open(input2,"r")do file2
            k = 0
            for line in eachline(file)
                R2=readline(file2)
                if line[1] == '@' && mode == ""
                    header = R2
                    mode = "header"
                elseif mode == "header"
                    k+=1
                    #seq = line[21:length(line)]
                    seq = line
                    mode = "line3"
                    maximum_bbc_score=0
                    maximum_bbc_score_number=0
                    for (i, row) in enumerate( eachrow( bbc_df ) ) 
                        bbc_score = BBCSemiglobalAlignmentScore(seq,row.Full_seq,maximum_errorrate)
                        if maximum_bbc_score<bbc_score
                            maximum_bbc_score = bbc_score
                            maximum_bbc_score_number = i
                        end
                    end
                    seq = R2
                    if maximum_bbc_score_number == 0
                        filename = output_dir * "/trimmed-unknown.fastq"
                    else 
                        filename = output_dir*"/trimmed-"*string(bbc_df.ID[maximum_bbc_score_number])*".fastq"    
                    end
                elseif mode == "line3"
                    line3 = R2
                    mode = "quality_score"
                elseif mode =="quality_score"
                    quality_score = R2
                    open(filename,"a")do outputfile
                        write(outputfile,header*"\n"*seq*"\n"*line3*"\n"*quality_score*"\n")
                    end 
                    mode = ""
                end
            end
        end
    end
end
end

