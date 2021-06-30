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

   # DOS assingments
   prm = AnnA.load_parameters(;N=200)
   prm.nᵢ=1e11u"cm^-3"
   prm.gcₑ= prm.gc 
   prm.gvₕ= prm.gv 
   @test prm.gc  ≈ prm.gv  ≈ prm.gcₑ  ≈ prm.gvₕ ≈ 2.750508782111764e30u"m^-3"
   
end