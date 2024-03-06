function recoMP2RAGE_CS(b::BrukerFile;λ=0.01,iteration=30, zeroFilled = false)

    raw = RawAcquisitionData_MP2RAGE_CS(b);
    acq = AcquisitionData(raw,OffsetBruker = true);

    if b["EnableCS"] == "Fully" || zeroFilled
        @info "FFT reconstruction"
        params = Dict{Symbol, Any}()
        params[:reco] = "direct"
        params[:reconSize] = tuple(acq.encodingSize...)
        Ireco = reconstruction(acq, params);

        sens = estimateCoilSensitivities(Ireco[:,:,:,2:2,:], 1.e-2)
        Icomb = sum(conj.(sens).*Ireco,dims=5) ./sum(abs.(sens).^2,dims = 5)

    else
        @info "CS reconstruction"
        kspace = extract3DKSpace(acq);

        calibSize = parse.(Int,b["CenterMaskSize"])
        calibData = crop(kspace[:,:,:,2,:],(calibSize,calibSize,calibSize))
        kspace = Nothing # clear data
        sensSize = (acq.encodingSize[1],acq.encodingSize[2],acq.encodingSize[3]) # (size(kspace,1),size(kspace,2),size(kspace,3))
        sens = espirit(calibData,sensSize ,(6,6,6),eigThresh_2 = 0); # really long on my mac

        sens = reshape(sens,acq.encodingSize[1],acq.encodingSize[2],acq.encodingSize[3],:);

        # CS reco
        T = eltype(real(sens[1]))
        params = Dict{Symbol, Any}()
        params[:reco] = "multiCoil"
        params[:reconSize] = sensSize
        params[:senseMaps] = sens;

        params[:solver] = "admm"
        params[:sparseTrafoName] = "Wavelet"
        params[:regularization] = "L1"
        params[:λ] = T(λ) # 5.e-2
        params[:iterations] = iteration
        params[:normalize_ρ] = true
        params[:ρ] = T(0.95)
        params[:relTol] = 0.1
        params[:normalizeReg] = true

        I_wav = reconstruction(acq, params);
        Icomb = I_wav.data
    end

    # apply MP2RAGE
    MP2 = mp2rage_comb(Icomb[:,:,:,:,1])

    TI₁ = parse(Float64,b["EffectiveTE"][1])
    TI₂ = parse(Float64,b["EffectiveTE"][2])
    TR = parse(Float64,b["PVM_RepetitionTime"])
    MP2RAGE_TR = parse(Float64,b["RecoveryTime"])
    ETL = parse(Int64,b["RareFactor"])
    α₁ = parse.(Float64,split(b["ExcPulse1"],", ")[3])
    α₂ = parse.(Float64,split(b["ExcPulse2"],", ")[3])

    params_MP2 = ParamsMP2RAGE(TI₁,TI₂,TR,MP2RAGE_TR,ETL,α₁,α₂)
    T1map, = mp2rage_T1maps(MP2,params_MP2)

    return Icomb,MP2,T1map,sens,acq
end


"""
Documentation for RawAcquisitionData_MP2RAGE_CS(b::BrukerFile)
"""
function RawAcquisitionData_MP2RAGE_CS(b::BrukerFile)
    # extract parameters
    sx = MRIFiles.acqSize(b)[1]
    sy = MRIFiles.pvmMatrix(b)[2]
    sz = MRIFiles.pvmMatrix(b)[3]
    RareFact = parse.(Int,b["RareFactor"])

    # extract position of rawdata
    GradPhaseVect = parse.(Float32,b["GradPhaseVector"])
    GradPhaseVect=GradPhaseVect*sy/2;
    GradPhaseVect=round.(Int,(GradPhaseVect.-minimum(GradPhaseVect)));

    GradSliceVect = parse.(Float32,b["GradSliceVector"])
    GradSliceVect=GradSliceVect*sz/2;
    GradSliceVect=round.(Int,(GradSliceVect.-minimum(GradSliceVect)))

    #=
    mask = zeros(Int,(sy,sz))
    anim = @animate for i = 1 : 100#length(GradPhaseVect)
        mask[GradPhaseVect[i],GradSliceVect[i]]=1;
        heatmap(mask,color = :greys,aspect_ratio = 1,legend = :none , axis=nothing,lims=(1,sz))
    end
    gif(anim, "anim_fps15.gif", fps = 15)
    =#

    # read data
    dtype = Complex{MRIFiles.acqDataType(b)}
    filename = joinpath(b.path, "fid")

    N = MRIFiles.acqSize(b)
    numChannel = parse.(Int,b["PVM_EncNReceivers"])
    profileLength = Int((ceil(N[1]*numChannel*sizeof(dtype)/1024))*1024/sizeof(dtype)) # number of points + zeros

    RareFact = MRIFiles.acqRareFactor(b)
    numSlices = MRIFiles.acqNumSlices(b)
    numEchos = MRIFiles.acqNumEchos(b)
    numEncSteps2 = length(N) == 3 ? N[3] : 1
    numRep = MRIFiles.acqNumRepetitions(b)

    I = open(filename,"r") do fd
        read!(fd,Array{dtype,5}(undef, profileLength,
                                        RareFact,
                                       numEchos,
                                       div(N[2], RareFact),
                                       numRep))[1:N[1]*numChannel,:,:,:,:]
      end;

    I2 = permutedims(I,(1,2,4,3,5)); # put echos in last dimension
    I2 = reshape(I2,sx,numChannel,N[2],numEchos,numRep)


    objOrd = MRIFiles.acqObjOrder(b)
    objOrd = objOrd.-minimum(objOrd)

    gradMatrix = MRIFiles.acqGradMatrix(b)

    offset1 = MRIFiles.acqReadOffset(b)
    offset2 = MRIFiles.acqPhase1Offset(b)
    offset3 = ndims(b) == 2 ? MRIFiles.acqSliceOffset(b) : MRIFiles.pvmEffPhase2Offset(b)

    profiles = MRIBase.Profile[]
    for nR = 1:numRep
        for nEcho=1:numEchos
            for nEnc = 1:N[2]

            counter = EncodingCounters(kspace_encode_step_1=GradPhaseVect[nEnc],
            kspace_encode_step_2=GradSliceVect[nEnc],
            average=0,
            slice=0,
            contrast=nEcho-1,
            phase=0,
            repetition=nR-1,
            set=0,
            segment=0)

            nSl = 1
            G = gradMatrix[:,:,nSl]
                    read_dir = (G[1,1],G[2,1],G[3,1])
                    phase_dir = (G[1,2],G[2,2],G[3,2])
                    slice_dir = (G[1,3],G[2,3],G[3,3])

            # Not sure if the following is correct...
            pos = offset1[nSl]*G[:,1] +
                    offset2[nSl]*G[:,2] +
                    offset3[nSl]*G[:,3]

            position = (pos[1], pos[2], pos[3])

            head = AcquisitionHeader(number_of_samples=sx, idx=counter,
                                        read_dir=read_dir, phase_dir=phase_dir,
                                        slice_dir=slice_dir, position=position,
                                        center_sample=div(sx,2),
                                        available_channels = numChannel, #TODO
                                        active_channels = numChannel)
            traj = Matrix{Float32}(undef,0,0)
            dat = map(ComplexF32, reshape(I2[:,:,nEnc,nEcho,nR],:,numChannel))
            push!(profiles, MRIBase.Profile(head,traj,dat) )
            end
        end
    end
    params = Dict{String,Any}()
    params["trajectory"] = "Cartesian"
    params["encodedSize"] = [sx;sy;sz]
    F = MRIFiles.acqFov(b)
    params["encodedFOV"] = F*1000.0 #convert m to mm
    params["receiverChannels"] = numChannel
    params["H1resonanceFrequency_Hz"] = round(Int, parse(Float64,b["SW"])*1000000)
    params["studyID"] = b["VisuStudyId"]
    #params["studyDescription"] = b["ACQ_scan_name"]
    #params["studyInstanceUID"] =
    params["referringPhysicianName"] = MRIFiles.latin1toutf8(b["ACQ_operator"])

    params["patientName"] = b["VisuSubjectName"]

    params["measurementID"] = parse(Int64,b["VisuExperimentNumber"])
    params["seriesDescription"] = b["ACQ_scan_name"]

    params["institutionName"] = MRIFiles.latin1toutf8(b["ACQ_institution"])
    params["stationName"] = b["ACQ_station"]
    params["systemVendor"] = "Bruker"

    params["TR"] = parse.(Float32,b["RecoveryTime"])
    params["delta_TR"] = MRIFiles.acqRepetitionTime(b)
    params["TE"] = parse.(Float32,b["PVM_EchoTime"])
    params["TI"] = parse.(Float32,b["EffectiveTE"])
    params["flipAngle_deg"] = parse.(Int64,[split(b["ExcPulse1"], ", ")[3];split(b["ExcPulse2"], ", ")[3]])
    params["sequence_type"] = MRIFiles.acqProtocolName(b)
    params["echo_spacing"] = MRIFiles.acqInterEchoTime(b)

    raw = RawAcquisitionData(params, profiles);
    return raw
end

"""
Documentation extractKSpace
"""
function extract3DKSpace(acqData::AcquisitionData{T}) where T
    if !MRIReco.isCartesian(trajectory(acqData, 1))
        @error "espirit does not yet support non-cartesian sampling"
    end
    nx, ny, nz = acqData.encodingSize[1:3]
    numChan, numSl = MRIFiles.numChannels(acqData), MRIFiles.numSlices(acqData)
    numEcho = length(acqData.traj)
    kdata = zeros(Complex{T}, nx * ny * nz, numEcho,numChan)

    for echo = 1:numEcho
        for coil = 1:numChan
            kdata[acqData.subsampleIndices[numEcho], echo,coil] .= kData(acqData, echo, coil, 1)
        end
    end
    kdata = reshape(kdata, nx, ny, nz, numEcho, numChan)
    return kdata
end

"""
Documentation : Crop the central area
"""
function crop(A::Array{T,4}, s::NTuple{3,Int64}) where {T}
    nx, ny, nz = size(A)
    idx_x = div(nx, 2)-div(s[1], 2)+1:div(nx, 2)-div(s[1], 2)+s[1]
    idx_y = div(ny, 2)-div(s[2], 2)+1:div(ny, 2)-div(s[2], 2)+s[2]
    idx_z = div(nz, 2)-div(s[3], 2)+1:div(nz, 2)-div(s[3], 2)+s[3]
    return A[idx_x, idx_y, idx_z,:]
end
