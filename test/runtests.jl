using DemultiplexJulia
using Test

@testset "DemultiplexJulia.jl" begin
    # Write your tests here.
    @test BBCSemiglobalAlignmentScore("AAAAAGGAA","GG",0.22) == 2
    @test BBCSemiglobalAlignmentScore("TTTTGGGGCCCGAAATGTGAAATAAGGGTA","CCCGAA",0.22) == 6
    @test BBCSemiglobalAlignmentScore("TATATTCGCGGCATCCGCACTTATGCGCGAATGCATAGCCAACTGCTTGGCAGTTGGCTTCAGACCA","AGCCAACTGCTTGGCAGTTGGC",0.22) == 22
    @test BBCSemiglobalAlignmentScore("TATATTCGCGGCATCCGCACTTATGCGCGAATGCATAGCCAACTGCTTGGCAGTTGGCTTCAGACCA","AGCCAACTGCTTGGCAGTTGGC",0.22) == 22    
end
