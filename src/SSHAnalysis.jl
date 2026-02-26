module SSHAnalysis
using LinearAlgebra
using QuadGK
using DelimitedFiles
using Plots
using LaTeXStrings
using ProgressMeter
using Base.Threads
using Dates

export GenerateDataEntropy, GenerateDataEntropyEquilibrium, 
       GenerateDataVarianceEquilibrium, GenerateDataVariance,  
       SaveWithTimeStamp, ProjectRoot, RealSpaceDensityM,
       @showprogress, EntanglementEntropy, SingleVariance,
       EntropyandVariance_sub, EigenvaluesDensity_sub

include("BlochVectors.jl")
include("DensityMatrix.jl")
include("Entropy.jl")
include("Variance.jl")
include("Utils.jl")

end
