using DemultiplexJulia
using Test


function test_Demultiplex()
    # Test inputs
    tmpdir = mktempdir()

    input = joinpath(tmpdir, "test_input.fastq")
    input2 = joinpath(tmpdir, "test_input2.fastq")
    bbc_tsv = joinpath(tmpdir, "test_bbc.tsv")
    output_dir = joinpath(tmpdir, "test_output")

    open(input, "w") do io
        println(io, "@example_read1")
        println(io, "AATTATATGGGGGGATATT")
        println(io, "+")
        println(io, rand('A':'J', 19))
        println(io, "@example_read2")
        println(io, "AATTATATGCCGCCATATT")
        println(io, "+")
        println(io, rand('A':'J', 19))
    end

    open(input2, "w") do io
        println(io, "@example_read1")
        println(io, "ATGC")
        println(io, "+")
        println(io, "FFFF")
        println(io, "@example_read2")
        println(io, "TTTTT")
        println(io, "+")
        println(io, "FFF:F")
    end

    io = open(bbc_tsv, "w")
    data = [
        ["ID", "Number", "Full_seq", "Full_annotation"],
        ["sample_1", 1, "AAAGGGGGGAAAT", "LLLBBBBBB3333"],
        ["sample_2", 2, "AACCCCCCAAA", "LLBBBBBB333"],
        ["sample_3", 3, "AAAGGCGGCAT", "LLLBBBBBB33"]
    ]
    for row in data
        row_str = join(string.(row), "\t")
        println(io, row_str)
    end
    close(io)

    # Run function
    Demultiplex(input, input2, bbc_tsv, output_dir, 0.22)

    # Test outputsq
    @test isdir(output_dir)
    @test !isfile(joinpath(output_dir, "trimmed-unknown.fastq"))
    @test !isfile(joinpath(output_dir, "trimmed-sample_1.fastq"))
    @test isfile(joinpath(output_dir, "trimmed-sample_2.fastq"))
    @test isfile(joinpath(output_dir, "trimmed-sample_3.fastq"))
    open(joinpath(output_dir, "trimmed-sample_2.fastq"),"r") do out
        @test readline(out) == "@example_read1"
        @test readline(out) == "ATGC"
        @test readline(out) == "+"
        @test readline(out) == "FFFF"
    end
    open(joinpath(output_dir, "trimmed-sample_3.fastq"),"r") do out
        @test readline(out) == "@example_read2"
        @test readline(out) == "TTTTT"
        @test readline(out) == "+"
        @test readline(out) == "FFF:F"
    end
end

function test_BBCSemiglobalAlignmentScore()
    ref = "AAAAAGGAA"
    query = "GG"
    maximum_errorrate = 0.22

    score = BBCSemiglobalAlignmentScore(ref, query, maximum_errorrate)
    @test score == 2

    ref = "TTTTGGGGCCCGAAATGTGAAATAAGGGTA"
    query = "CCCGAA"
    score = BBCSemiglobalAlignmentScore(ref, query, maximum_errorrate)
    @test score == 6

    ref = "TATATTCGCGGCATCCGCACTTATGCGCGAATGCATAGCCAACTGCTTGGCAGTTGGCTTCAGACCA"
    query = "AGCCAACTGCTTGGCAGTTGGC"
    score = BBCSemiglobalAlignmentScore(ref, query, maximum_errorrate)
    @test score == 22

    ref = "TATATTCGCGGCATCCGCACTTATGCGCGAATGCATAGCCAACTGCTTGCAGTTGGCTTCAGACCA"
    score = BBCSemiglobalAlignmentScore(ref, query, maximum_errorrate)
    @test score == 19
end

@testset "DemultiplexJulia.jl" begin
    @testset "Demultiplex tests" begin
        test_Demultiplex()
    end

    @testset "BBCSemiglobalAlignmentScore tests" begin
        test_BBCSemiglobalAlignmentScore()
    end
end
