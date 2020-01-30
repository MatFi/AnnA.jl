mutable struct OCVDProblem{B<:Bool,T<:Number,P<:Parameters,C<:AlgControl} <: AbstactProblem
    parameters::P
    on_time::T
    decay_time::T
    alg_control::C
end
