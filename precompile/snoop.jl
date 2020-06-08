using SnoopCompile

sol_problems = @snoopi tmin=0.01 include("./scripts.jl")
pc = SnoopCompile.parcel(sol_problems)
SnoopCompile.write("./pc", pc)
