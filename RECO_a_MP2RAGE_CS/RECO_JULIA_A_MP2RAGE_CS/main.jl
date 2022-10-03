using Pkg
Pkg.instantiate()

using MRIReco
using qMRI
include("utils_a_MP2RAGE_CS.jl")

b = BrukerFile("data/20220930_144024_AT_MP2RAGE_JULIA_05_1_1/5")

Ireco,MP2,T1map,acq,sens = recoMP2RAGE_CS(b, λ=0.1,iteration=30,zeroFilled = false)

# write NIFTI
T1map = reshape(T1map,collect(size(T1map))...,1,1)
T1map = MRIReco.makeAxisArray(Float64.(T1map),acq)
saveImage("data/test.nii", T1map)
