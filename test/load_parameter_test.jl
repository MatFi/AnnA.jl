@testset "load_parameters" begin
   @test AnnA.load_parameters(;N=200) isa AbstractParameters
   
   #test for failure
   cpath= pwd()
   cd("..")
   @test_throws LoadError AnnA.load_parameters(;N=200)
   write_template()
   @test AnnA.load_parameters(;N=200) isa AbstractParameters
   rm("Parameters.jl")
   cd(cpath)
end