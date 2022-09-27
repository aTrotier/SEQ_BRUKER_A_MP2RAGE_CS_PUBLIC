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
            lines2[i]="##\$VisuCoreDataMin= ( $nframes )"
            str = prod([string(VisuCoreDataMin[i]) * " " for i=1:nframes])
            str = str[1:end-1]
            lines2[i+1]="$str"
        end

        if contains(lines[i],"##\$VisuCoreDataMax=")
            lines2[i]="##\$VisuCoreDataMax= ( $nframes )"
            str = prod([string(VisuCoreDataMax[i]) * " " for i=1:nframes])
            str = str[1:end-1]
            lines2[i+1]="$str"
        end

        if contains(lines[i],"##\$VisuCoreDataOffs=")
            lines2[i]="##\$VisuCoreDataOffs= ( $nframes )"
            str = prod([string(VisuCoreDataOffs[i]) * " " for i=1:nframes])
            str = str[1:end-1]
            lines2[i+1]="$str"
        end

        if contains(lines[i],"##\$VisuCoreDataSlope=")
            lines2[i]="##\$VisuCoreDataSlope= ( $nframes )"
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
