%% Script : Compare CS reconstruction to groundTruth
%% Get Bruker Acquistion parameters
s_in = OBJ_MP2RAGE_RECO; % get the path to bruker file
bruk = read_bru_experiment(s_in.BRUKER_PATH); % read parameter

%% Create Mask
% Set MASK parameters
% [mask] = poissonDiscMex(fovy,fovz,sky,skz,ry,rz,ncal,cutcorners, pp)

% parameter fixed by the acquisition
fovy = bruk.method.PVM_Fov(2);
fovz = bruk.method.PVM_Fov(3);
sky = bruk.method.PVM_Matrix(2);
skz = bruk.method.PVM_Matrix(3);

% Free parameter to create the mask
bruk.method.CenterMaskSize = 24; % TO REMOVE
ncal = bruk.method.CenterMaskSize; % pour le moment on garde la meme valeur qu'en fully
ry = 2;
rz = 2;
cutcorners = 0;
pp = 1;

% Create mask
% If not present -> create the c library with command : mex poissonDiscMex.c
[mask] = poissonDiscMex(fovy,fovz,sky,skz,ry,rz,ncal,cutcorners, pp);
figure;imagesc(mask)

% Number of point in the mask
maskPoint = sum(mask(:))

Accel = length(mask(:))/maskPoint

%% Create fully and undersample dataset :

% create kspace
kdata_bart = Bruker_to_buffer(bruk);

if strcmp(s_in.OPTION_RECO,'K_AVERAGE')
    % kdata_bart = mean(kdata_bart,8);
    tmp = sum(kdata_bart~=0,8);
    tmp(tmp == 0) = 1;
    kdata_bart = sum(kdata_bart,8) ./ tmp;
end


% Undersampled dataset : 

mask_kdata = repmat(permute(mask,[3 1 2]),[size(kdata_bart,1) 1 1 size(kdata_bart,4) 1 2]);
size(mask_kdata)

kdata_under = kdata_bart .* mask_kdata;




%% Find best reco for W in spatial direction and TV  
% 1) spatial dimension (CS)
% 2) TI dimension (CS2)

 reco_in.kdata_fully = kdata_bart;
 reco_in.kdata_under = kdata_under;
 
 reco_in.Wlist = [0 0.0001 0.001 0.005 0.0075 0.01:0.01:0.1 0.2:0.1:1];
 reco_in.TVlist = [0 0.0001 0.001 0.005 0.0075 0.01:0.01:0.1 0.2:0.1:1];
 reco_in.ITERlist = [200];

 reco_in.GPU = 1;
 reco_in.calib = ncal;
 reco_in.SENS_MAPS = 2;
 reco_in.bg_mult = 3;
 
 [s_out] = reco_multi_param(reco_in);
 
% What is the best reco :

[a,b]=min(s_out.RMSD(:));
[idx1,idx2,idx3] = ind2sub([size(s_out.RMSD,1) size(s_out.RMSD,2) size(s_out.RMSD,3)],b);
fprintf('RMSD = %f - TV = %f - W = %f - IT = %f \n',a, reco_in.TVlist(idx1),reco_in.Wlist(idx2),reco_in.ITERlist(idx3));


[a,b]=min(s_out.RMSD2(:));
[idx1,idx2,idx3] = ind2sub([size(s_out.RMSD2,1) size(s_out.RMSD2,2) size(s_out.RMSD2,3)],b);
fprintf('RMSD2 = %f - TV = %f - W = %f - IT = %f \n',a, reco_in.TVlist(idx1),reco_in.Wlist(idx2),reco_in.ITERlist(idx3));


[a,b]=max(s_out.SSIM(:));
[idx1,idx2,idx3] = ind2sub([size(s_out.SSIM,1) size(s_out.SSIM,2) size(s_out.SSIM,3)],b);
fprintf('SSIM = %f - TV = %f - W = %f - IT = %f \n',a, reco_in.TVlist(idx1),reco_in.Wlist(idx2),reco_in.ITERlist(idx3));

[a,b]=max(s_out.SSIM2(:));
[idx1,idx2,idx3] = ind2sub([size(s_out.SSIM2,1) size(s_out.SSIM2,2) size(s_out.SSIM2,3)],b);
fprintf('SSIM2 = %f - TV = %f - W = %f - IT = %f \n',a, reco_in.TVlist(idx1),reco_in.Wlist(idx2),reco_in.ITERlist(idx3));

% show best results
%imagine(s_out.MP2RAGE_fft,s_out.MP2RAGE_zf,s_out.MP2RAGE_CS,s_out.MP2RAGE_CS2,s_out.MP2RAGE_CS_SSIM,s_out.MP2RAGE_CS_SSIM2);

%% Find best reco for W in spatial direction only

 reco_inW.kdata_fully = kdata_bart;
 reco_inW.kdata_under = kdata_under;
 
 reco_inW.Wlist = [0 0.0001 0.001 0.005 0.0075 0.01:0.01:0.1 0.2:0.1:1];
 reco_inW.ITERlist = [200];

 reco_inW.GPU = 1;
 reco_inW.calib = ncal;
 reco_inW.SENS_MAPS = 2;
 
 [s_out_W] = reco_multi_param_W(reco_inW);
 
% What is the best reco for reco W only

[a,b]=min(s_out_W.RMSD(:))
[idx1,idx2,idx3] = ind2sub([size(s_out_W.RMSD,1) size(s_out_W.RMSD,2) size(s_out_W.RMSD,3)],b)
fprintf('RMSDE = %f - TV = %d - W = %d - IT = %d \n',a, reco_in.TVlist(idx1),reco_in.Wlist(idx2),reco_in.ITERlist(idx3));

% show best results
%imagine(s_out_W.im_fft,s_out_W.im_zf,s_out_W.im_CS,s_out_W.diff_CS,s_out.im_CS2,s_out.diff_CS2);

%% Do your own reco : 

% Set Reco Parameter
s_in.GPU = 1;
s_in.OPTION_RECO = "K_AVERAGE";
s_in.ITER = 200;
s_in.SENS_MAPS = 2;
s_in.W_LAMBDA = 0.03;
s_in.TV_LAMBDA = 0.1;


% sensitivity maps : check the results if needed add -c option to change the mask
if s_in.GPU
    sens_bart = bart(sprintf('ecalib -g -d3 -r %d -m %d',bruk.method.CenterMaskSize,s_in.SENS_MAPS), kdata_bart(:,:,:,:,:,2,1)); % second TI and first rep is used
else
    sens_bart = bart(sprintf('ecalib -d3 -r %d -m %d',bruk.method.CenterMaskSize,s_in.SENS_MAPS), kdata_bart(:,:,:,:,:,2,1)); % second TI and first rep is used
end
% Fully reco

im_fft = bart('pics -d4 -i1 -e',kdata_bart,sens_bart);
im_fft = abs(squeeze(im_fft(:,:,:,:,1,2)));
% Zero filling

im_zf = bart('pics -d4 -i1 -e',kdata_under,sens_bart);
im_zf = abs(squeeze(im_zf(:,:,:,:,1,2)));

% CS RECO W and TV in 3D

if s_in.GPU
    warning('GPU is used : if reconstruction failder try with GPU = 0');
    
   imCS = bart(sprintf('pics -d2 -g -U -i%d -e -R W:7:0:%d -R T:7:0:%d',s_in.ITER,s_in.W_LAMBDA,s_in.TV_LAMBDA),  kdata_under , sens_bart);
else
    imCS = bart(sprintf('pics -d2 -U -i%d -e -R W:7:0:%d -R T:7:0:%d',s_in.ITER,s_in.W_LAMBDA,s_in.TV_LAMBDA),  kdata_under , sens_bart);
end

if strcmp(s_in.OPTION_RECO,'IM_AVERAGE')
    imCS=squeeze(mean(imCS(:,:,:,:,1,:,:,:),8)); % use only first set of map
end

imCS = abs(squeeze(imCS(:,:,:,:,1,2)));

% CS RECO W 3D and TV in TI dim
if s_in.GPU
    warning('GPU is used : if reconstruction failder try with GPU = 0');
    
   imCS2 = bart(sprintf('pics -d2 -g -U -i%d -e -R W:7:0:%d -R T:32:0:%d',s_in.ITER,s_in.W_LAMBDA,s_in.TV_LAMBDA),  kdata_under , sens_bart);
else
    imCS2 = bart(sprintf('pics -d2 -U -i%d -e -R W:7:0:%d -R T:32:0:%d',s_in.ITER,s_in.W_LAMBDA,s_in.TV_LAMBDA),  kdata_under , sens_bart);
end

if strcmp(s_in.OPTION_RECO,'IM_AVERAGE')
    imCS2=squeeze(mean(imCS2(:,:,:,:,1,:,:,:),8)); % use only first set of map
end

imCS2 = abs(squeeze(imCS2(:,:,:,:,1,2)));

% Comparaison
diff_zf = (im_zf-im_fft);
diff_CS = (imCS-im_fft);
diff_CS2 = (imCS2-im_fft);

rms(diff_zf(:))
rms(diff_CS(:))
rms(diff_CS2(:))

%imagine(im_fft,im_zf,imCS,imCS2,diff_zf,diff_CS,diff_CS2);