using MRIReco

const BrukerWordType = ["_32BIT_SGN_INT","_16BIT_SGN_INT","_32BIT_FLOAT","_8BIT_UNSGN_INT"]


function my_bruker_write_in_REAL(image_4D::Array{Float32,4},dirname::String;calcOffsetSlope=true) 
    VisuCoreWordType = "_32BIT_FLOAT";
    VisuCoreByteOrder = "littleEndian";
    VisuCoreDataSlope = ones(size(image_4D,4));
    VisuCoreDataOffs = zeros(size(image_4D,4));
    VisuCoreFrameType ="REAL_IMAGE";

    VisuCoreSize=[size(image_4D,1), size(image_4D,2), size(image_4D,3), size(image_4D,4)];


    if(VisuCoreWordType == "_32BIT_SGN_INT")
        format="int32";
        round_type="round";
        MaxNorm=2^31-1;
        allowed_min=-2147483648 +1;
        allowed_max= 2147483647 -1;
    elseif (VisuCoreWordType == "_16BIT_SGN_INT")
        format="int16";
        round_type="round";
        MaxNorm=2^15-1;
        allowed_min=-32768 +1;
        allowed_max=32767 -1;
    elseif (VisuCoreWordType == "_32BIT_FLOAT")
        format="float32";
        round_type="single";
        allowed_min=-prevfloat(Inf32) +1;
        allowed_max=prevfloat(Inf32) -1;
    elseif (VisuCoreWordType == "_8BIT_UNSGN_INT")
        format="uint8";
        round_type="round";
        allowed_min=0 +1;
        allowed_max=255 -1;
    else
        @error "Data-Format not correct specified! VisuCoreWordType should be equal to : \n
        [\"_32BIT_SGN_INT\",\"_16BIT_SGN_INT\",\"_32BIT_FLOAT\",\"_8BIT_UNSGN_INT\"]" 
    end

    nframes = size(image_4D,4)

    # check if data is out of range
    VisuCoreDataMin = Float32[]
    VisuCoreDataMax = Float32[]
    for i = 1:nframes
	minVal = minimum(filter(!isnan,image_4D[:,:,:,i]))
	maxVal = maximum(filter(!isnan,image_4D[:,:,:,i]))
    
        (minVal < allowed_min) ? println("Some value of the image are inferior to the minimum allowed value") : 1
        (maxVal > allowed_max) ? println("Some value of the image are superior to the minimum allowed value") : 1
    
        image_4D = Float32.(image_4D)
        push!(VisuCoreDataMin,minVal)
        push!(VisuCoreDataMax,maxVal)
    end

    twodseq = open(dirname*"/2dseq","w");
    write(twodseq,image_4D)
    close(twodseq);

 
    f = open(dirname*"/visu_pars")
    lines = readlines(f)
    close(f)
    lines2 = copy(lines)

    for i = 1:length(lines)
        if contains(lines[i],"##\$VisuCoreFrameCount=") lines2[i]="##\$VisuCoreFrameCount=$nframes" end

        if contains(lines[i],"##\$VisuCoreDataMin=") 
            lines2[i]="##\$VisuCoreDataMin=( $nframes )" 
            str = prod([string(VisuCoreDataMin[i]) * " " for i=1:nframes])
            str = str[1:end-1]
            lines2[i+1]="$str"
        end

        if contains(lines[i],"##\$VisuCoreDataMax=") 
            lines2[i]="##\$VisuCoreDataMax=( $nframes )" 
            str = prod([string(VisuCoreDataMax[i]) * " " for i=1:nframes])
            str = str[1:end-1]
            lines2[i+1]="$str"
        end

        if contains(lines[i],"##\$VisuCoreDataOffs=") 
            lines2[i]="##\$VisuCoreDataOffs=( $nframes )" 
            str = prod([string(VisuCoreDataOffs[i]) * " " for i=1:nframes])
            str = str[1:end-1]
            lines2[i+1]="$str"
        end

        if contains(lines[i],"##\$VisuCoreDataSlope=") 
            lines2[i]="##\$VisuCoreDataSlope=( $nframes )" 
            str = prod([string(VisuCoreDataSlope[i]) * " " for i=1:nframes])
            str = str[1:end-1]
            lines2[i+1]="$str"
        end

        if contains(lines[i],"##\$VisuCoreFrameType=") lines2[i+1]="$VisuCoreFrameType" end

        if contains(lines[i],"##\$VisuCoreWordType=") lines2[i]="##\$VisuCoreWordType=$VisuCoreWordType" end

        if contains(lines[i],"##\$VisuFGOrderDesc=") 
            lines2[i+1]="($nframes, <FG_CYCLE>, <>, 0, 0)"
        end

        if contains(lines[i],"##\$VisuGroupDepVals=") 
            lines2[i+1]="(<>, 0)"
        end
    end

    path2visu = dirname*"/visu_pars"

    rm(path2visu,force=true)
    f = open(path2visu,"w")
    for i = 1:length(lines2)
        println(f,lines2[i])
    end
    close(f)
end





"""
        [VisuCoreDataSlope, VisuCoreDataOffs, VisuCoreDataMin, VisuCoreDataMax]=bruker_genSlopeOffsMinMax(VisuCoreWordType, data, VisuCoreDataSlope, VisuCoreDataOffs)

        Implementation in Julia of the function from bruker pv_matlab
            # TODO : 
            1) vérifier complex et matrix > 5 dims

% generates some important visu-parameters from a dataset, dependends on
% the target word-type.
% There are 2 modi:
%   - insert VisuCoreWordType, data, VisuCoreDataSlope=1 and VisuCoreDataOffs=0 -> will calculate all output
%     variables
%   - insert a different Slope or Offset -> will calculate only the VisuCoreDataMin,
%     VisuCoreDataMax based on data, VisuCoreDataSlope and VisuCoreDataOffs
%
% Note: VisuCoreWordType _32BIT_SGN_INT _16BIT_SGN_INT and _32BIT_FLOAT
%       and if you calculate all parameters, it will leave the zero-value by the zero-value, that means: 
%       it's possible that not the complete range is used.  
%
% IN: some Parameter of visu_pars, and the data-matrix with (dim1, dim2, dim3, dim4, NumberFrames)
% Out: some visu_pars parameter which are necessary for writing a PV-dataset 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright (c) 2012
% Bruker BioSpin MRI GmbH
% D-76275 Ettlingen, Germany
%
% All Rights Reserved
%
%  bruker_genSlopeOffsMinMax.m,v 1.1 2012/09/11 14:22:10 pfre Exp 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""
function bruker_genSlopeOffsMinMax(VisuCoreWordType::String, data::Array, VisuCoreDataSlope=[1], VisuCoreDataOffs=[0]) 

if(VisuCoreWordType == "_32BIT_SGN_INT")
    allowed_min=-2147483648 +1;
    allowed_max= 2147483647 -1;
    keepimagezero=true;
    calcOffsetSlope=true;
elseif (VisuCoreWordType == "_16BIT_SGN_INT")
    allowed_min=-32768 +1;
    allowed_max=32767 -1;
    keepimagezero=true;
    calcOffsetSlope=true;
elseif (VisuCoreWordType == "_32BIT_FLOAT")
    allowed_min=-prevfloat(Inf32) +1;
    allowed_max=prevfloat(Inf32) -1;
    keepimagezero=true;
    calcOffsetSlope=false;
elseif (VisuCoreWordType == "_8BIT_UNSGN_INT")
    allowed_min=0 +1;
    allowed_max=255 -1;
    keepimagezero=false;
    calcOffsetSlope=true;
else
    @error "Data-Format not correct specified! VisuCoreWordType should be equal to : \n
    [\"_32BIT_SGN_INT\",\"_16BIT_SGN_INT\",\"_32BIT_FLOAT\",\"_8BIT_UNSGN_INT\"]" 
end

# Read Mins and Maxs
dims = size(dataB)

if length(dims)<=5 && isreal(dataB[1])
    maxpos=zeros(size(dataB,5),1)
    maxneg=zeros(size(dataB,5),1)

    for i=1:size(data,5)
        maxpos[i]=maximum(dataB[:,:,:,:,i])
        maxneg[i]=minimum(dataB[:,:,:,:,i])
    end

elseif isreal(dataB[1])
    @warn "dims > 5 is not tested"
    frameCount = prod(dims[5:end])
    maxpos=zeros(frameCount,1);
    maxneg=zeros(frameCount,1);
    dimlist=dims[5:end];

    # generate mapping-table1:
    dimmatrix_tmp=zeros(prod(dimlist),length(dimlist));
    for i1=1:prod(dimlist)
        modulus=mod(i1-1,dimlist[end]);
        reminder=fix((i1-1)/dimlist[end]);
        dimmatrix_tmp[i1,length(dimlist)]=modulus+1;
        for i2=length(dimlist)-1:-1:1
            modulus=mod(reminder,dimlist[i2]);
            reminder=fix(reminder/dimlist[i2]);
            dimmatrix_tmp[i1,i2]=modulus+1;
        end      
    end
    # convert to 10 Dims:
    dimmatrix=ones(size(dimmatrix_tmp,1),10);
    dimmatrix[:,1:size(dimmatrix_tmp,2)]=dimmatrix_tmp;
    dimmatrix_tmp = nothing

    # generate mapping-table2:
    map=reshape(1:frameCount, dimlist);
    
    for i=1:size(dimmatrix,1)
        counter=map[dimmatrix[i,1], dimmatrix[i,2], dimmatrix[i,3], dimmatrix[i,4], dimmatrix[i,5], dimmatrix[i,6], dimmatrix[i,7], dimmatrix[i,8], dimmatrix[i,9], dimmatrix[i,10]];
        maxpos[counter]=maximum(data[:,:,:,:,dimmatrix[i,1], dimmatrix[i,2], dimmatrix[i,3], dimmatrix[i,4], dimmatrix[i,5], dimmatrix[i,6], dimmatrix[i,7], dimmatrix[i,8], dimmatrix[i,9], dimmatrix[i,10]]);
        maxneg[counter]=minimum(data[:,:,:,:,dimmatrix[i,1], dimmatrix[i,2], dimmatrix[i,3], dimmatrix[i,4], dimmatrix[i,5], dimmatrix[i,6], dimmatrix[i,7], dimmatrix[i,8], dimmatrix[i,9], dimmatrix[i,10]]);
    end

elseif ~isreal(dataB[1])
    @warn "dims > 5 is not tested"
    frameCount=prod(dims[5:end]);
    maxpos=zeros(2*frameCount,1);
    maxneg=zeros(2*frameCount,1);
    dimlist=dims[5:end];

    # generate mapping-table1:
    dimmatrix_tmp=zeros(prod(dimlist),length(dimlist));
    for i1=1:prod(dimlist)
        modulus=mod(i1-1,dimlist[end]);
        reminder=fix((i1-1)/dimlist[end]);
        dimmatrix_tmp[i1,length(dimlist)]=modulus+1;
        for i2=length(dimlist)-1:-1:1
            modulus=mod(reminder,dimlist[i2]);
            reminder=fix(reminder/dimlist[i2]);
            dimmatrix_tmp[i1,i2]=modulus+1;
        end      
    end
    # convert to 10 Dims:
    dimmatrix=ones(size(dimmatrix_tmp,1),10);
    dimmatrix[:,1:size(dimmatrix_tmp,2)]=dimmatrix_tmp;
    dimmatrix_tmp = nothing

    # generate mapping-table2:
    map=reshape(1:frameCount, [dimlist,1]);
    
    for i=1:size(dimmatrix,1)
        counter=map[dimmatrix[i,1], dimmatrix[i,2], dimmatrix[i,3], dimmatrix[i,4], dimmatrix[i,5], dimmatrix[i,6], dimmatrix[i,7], dimmatrix[i,8], dimmatrix[i,9], dimmatrix[i,10]];
        maxpos[counter]=maximum(real(data[:,:,:,:,dimmatrix[i,1], dimmatrix[i,2], dimmatrix[i,3], dimmatrix[i,4], dimmatrix[i,5], dimmatrix[i,6], dimmatrix[i,7], dimmatrix[i,8], dimmatrix[i,9], dimmatrix[i,10]]));
        maxneg[counter]=minimum(real(data[:,:,:,:,dimmatrix[i,1], dimmatrix[i,2], dimmatrix[i,3], dimmatrix[i,4], dimmatrix[i,5], dimmatrix[i,6], dimmatrix[i,7], dimmatrix[i,8], dimmatrix[i,9], dimmatrix[i,10]]));
        maxpos[counter+frameCount]=maximum(imag(data[:,:,:,:,dimmatrix[i,1], dimmatrix[i,2], dimmatrix[i,3], dimmatrix[i,4], dimmatrix[i,5], dimmatrix[i,6], dimmatrix[i,7], dimmatrix[i,8], dimmatrix[i,9], dimmatrix[i,10]]));
        maxneg[counter+frameCount]=minimum(imag(data[:,:,:,:,dimmatrix[i,1], dimmatrix[i,2], dimmatrix[i,3], dimmatrix[i,4], dimmatrix[i,5], dimmatrix[i,6], dimmatrix[i,7], dimmatrix[i,8], dimmatrix[i,9], dimmatrix[i,10]]));
    end
end

## Claculate slope and off
# later in write2dseq: discdata=(data-offset)/slope

if calcOffsetSlope
    if prod(VisuCoreDataSlope[:])==1 && sum(abs.(VisuCoreDataOffs[:]))==0
        VisuCoreDataOffs=zeros(length(maxpos));
        VisuCoreDataSlope=zeros(length(maxpos));
        for i=1:length(maxpos)
            if keepimagezero
                slope_pos=maxpos[i]/allowed_max;
                slope_neg=maxneg[i]/allowed_min;
                VisuCoreDataSlope[i]=maximum([slope_pos, slope_neg]);
                # VisuCoreDataOffs=0; (still done)
            else
                VisuCoreDataSlope[i]=(maxpos[i]-maxneg[i])/(allowed_max-allowed_min);      
                VisuCoreDataOffs[i]=maxneg[i]-allowed_min*VisuCoreDataSlope[i];
            end
        end
    end
else
    VisuCoreDataOffs=zeros(length(maxpos),1);
    VisuCoreDataSlope=ones(length(maxpos),1);
end

## Calculate VisuCoreDataMax and VisuCoreDataMin
VisuCoreDataMax=zeros(length(maxpos));
VisuCoreDataMin=zeros(length(maxpos));
for i=1:length(maxpos)
    VisuCoreDataMax[i]=(maxpos[i]-VisuCoreDataOffs[i])/VisuCoreDataSlope[i];
    VisuCoreDataMin[i]=(maxneg[i]-VisuCoreDataOffs[i])/VisuCoreDataSlope[i];
    if VisuCoreDataMax[i] > allowed_max+1
        @warn "Your dataslope and data of are cutting the values"
        VisuCoreDataMax=allowed_max+1;
    end
    if VisuCoreDataMin[i] < allowed_min-1
        @warn "Your dataslope and data of are cutting the values"
    end
end

return VisuCoreDataSlope, VisuCoreDataOffs, VisuCoreDataMin, VisuCoreDataMax
end



