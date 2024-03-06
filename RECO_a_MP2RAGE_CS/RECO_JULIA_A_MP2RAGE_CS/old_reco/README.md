
## Online Reconstruction (NOT MAINTAINED anymore moved into old_reco folder)
The reconstruction can be performed either directly on the Bruker scanner by calling a macro. The results will be available in a 2dseq files and can be seen directly in Paravision.
Or you can reconstruct offline the data and write directly the nifti (without any information about the FOV position etc)

Two packages are required :
- `MRIReco` : read bruker files and perform advanced reconstruction
- `qMRI` : process the image and create T1 maps

Some functions are also requires and stores in 2 files :
- `utils_a_MP2RAGE_CS.jl` : Reorder the fid data according to the Compressed Sensing sampling patterns. Additionnal functions to extract the kspace and crop the center.
- `utils_write_bruker.jl` : write the resulting image in the 2dseq files and edit the `VISU_PARS`

## Launch Julia
Type `julia` in your terminal. The most important options you can add are :
- `-t auto` or `-t N` où N est le nombre de threads que vous voulez utiliser
- `--project PATH` let julia read the project and manifest in the directory **PATH**



## Offline reconstruction
In order to perform an offline reconstruction (without writing to the 2dseq folder) use the `main.jl` script.
Just edit the path to your experiment :
```julia
b = BrukerFile(PATH)
```

If you want to write the file to nifti you can run the 3 last lines.

