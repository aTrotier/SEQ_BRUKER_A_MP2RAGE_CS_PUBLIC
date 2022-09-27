using Pkg
#Pkg.activate(@__DIR__)
Pkg.instantiate()

using MRIReco
using qMRI
using Logging

disable_logging(LogLevel(10000))

#on bruker don't forget to create a file : ~/.julia/config/startup.jl
# ENV["HTTP_PROXY"] = "http://192.168.10.1:3128"
# ENV["HTTPS_PROXY"] = "http://192.168.10.1:3128"

println("Script Path = "*ARGS[1])
if isfile(ARGS[1]*"/2dseq")
  path_reco = ARGS[1]
  path_fid = ARGS[1]*"/../.."
else
  path_reco = ARGS[1]*"/pdata/1"
  path_fid = ARGS[1]
end
println("path_reco = "*path_reco)
println("path_fid = "*path_fid)

b=MRIReco.BrukerFile(path_fid)

# Reconstruction of a_MP2RAD
include("../utils_a_MP2RAGE_CS.jl");
Isense,MP2_sense,T1map = recoMP2RAGE_CS(b)


## Write back data in Bruker format
include("utils_write_bruker.jl");

image_4D = zeros(Float64,size(MP2_sense,1),size(MP2_sense,2),size(MP2_sense,3),4);
image_4D[:,:,:,1:2]=abs.(Isense[:,:,:,1:2])
image_4D[:,:,:,3]=MP2_sense;
image_4D[:,:,:,4]=T1map;

image_4D = Float32.(image_4D)

# Permute according to direction
@info b["PVM_SPackArrSliceOrient"][1]
@info b["PVM_SPackArrReadOrient"][1]

if b["PVM_SPackArrSliceOrient"][1] == "axial" && b["PVM_SPackArrReadOrient"][1] == "A_P"
  image_4D = permutedims(image_4D, (2,1,3,4))
elseif b["PVM_SPackArrSliceOrient"][1] == "coronal" && b["PVM_SPackArrReadOrient"][1] == "H_F"
  image_4D = permutedims(image_4D, (2,1,3,4))
elseif b["PVM_SPackArrSliceOrient"][1] == "sagittal" && b["PVM_SPackArrReadOrient"][1] == "H_F"
  image_4D = permutedims(image_4D, (2,1,3,4))
end


my_bruker_write_in_REAL(image_4D,path_reco)

exit(0)
