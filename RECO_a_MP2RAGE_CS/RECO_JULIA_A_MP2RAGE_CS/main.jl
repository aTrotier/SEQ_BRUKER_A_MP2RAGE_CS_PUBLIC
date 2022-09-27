using Pkg
Pkg.instantiate()

using MRIReco
using qMRI
include("utils_a_MP2RAGE_CS.jl")

b = BrukerFile("data/20220927_091210_AT_MP2RAGE_JULIA_04_1_1/3")
Ireco,MP2,T1map,acq = recoMP2RAGE_CS(b,forceCS = true)

# write NIFTI
T1map = reshape(T1map,collect(size(T1map))...,1,1)
T1map = MRIReco.makeAxisArray(Float64.(T1map),acq)
saveImage("data/test.nii", T1map)
