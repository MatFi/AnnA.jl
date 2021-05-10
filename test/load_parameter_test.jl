@testset "load_parameters" begin
   #should fail in test directory
   @test_throws LoadError load_parameters(;N=200)
   #but should wour if we create the template
   write_template()
   @test load_parameters(;N=200) isa AnnA.AbstractParameters
   rm("Parameters.jl")
   cpath= pwd()
   cd("..")
   # it should reflect the cd()
   @test AnnA.load_parameters(;N=200) isa AnnA.AbstractParameters

end