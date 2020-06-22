% imOut=a_MP2RAGE_CS(s)
%
% Author:   Aurélien TROTIER  (a.trotier@gmail.com)
% Date:     2020-04-16
% Partner: none
% Institute: CRMSB (Bordeaux, FRANCE)
%
% Function description:
%       Perform a CS reconstruction using BART (Berkley Advanced
%       Reconstruction ToolBox) for the MP2RAGE sequence acquired on Bruker
%       scanner.
%
% Input:
%       s_in : class generated from OBJ_MP2RAGE_RECO
%
% Output: 
%       s_out.T1map : T1 mapping in ms
%            .LUT : substructure for generation of lookup table
%            .MP2RAGE_mask : MP2RAGE with background removal (the MP2RAGE with background can be found under .LUT.MP2RAGE
%            .imTI : image of each invertion time
%            .param_in : OBJ_MP2RAGE_RECO used as input
%
% Algorithm & bibliography:
%        Trotier AJ, Rapacchi S, Faller TL, Sylvain M, Ribot EJ. 
%        Compressed-Sensing MP2RAGE sequence: application to the 
%        detection of brain metastases in mice at 7T. Magn Reson Med. 2018;(June):19. 
%
% See also :
%
% TO DO :
%   - Manage "K_AVERAGE if mask is different"

function [s_out]=a_MP2RAGE_CS_bart(s_in)
tic
if ~strcmp(class(s_in),"OBJ_MP2RAGE_RECO")
   error("Input is not a OBJ_MP2RAGE_RECO class object : Defined it using -> s_in = OBJ_MP2RAGE_RECO");
end
disp(s_in);
%% Read data and sequence parameters 
bruk = read_bru_experiment(s_in.BRUKER_PATH);

%%
sx = bruk.acqp.ACQ_size(1)/2;
sy = bruk.method.PVM_Matrix(2);
sz = bruk.method.PVM_Matrix(3);

GradPhaseVector = bruk.method.GradPhaseVector;
GradSliceVector = bruk.method.GradSliceVector;

GradPhaseVector=GradPhaseVector*sy/2;
GradPhaseVector=round(GradPhaseVector-min(GradPhaseVector)+1);
GradSliceVector=GradSliceVector*sz/2;
GradSliceVector=round(GradSliceVector-min(GradSliceVector)+1);

Nc=bruk.method.PVM_EncNReceivers;
NEchos=bruk.method.PVM_NEchoImages;
RareFactor=bruk.method.RareFactor;
NR=bruk.acqp.NR;

Genpoint = bruk.method.Genpoint;

%% mask
mask=zeros(sy,sz);

mask(sub2ind(size(mask),GradPhaseVector(:),GradSliceVector(:)))=1;
figure;imshow(mask);

%% read fid file in the same directory
rawData=bruk.fid;
% Remove 0 from data
sizeR2=128*ceil(sx*Nc/128);
rawData=reshape(rawData,sizeR2,[]);
rawData=rawData(1:Nc*sx,:);

rawDataTmp=cell(Nc,1);
for i=1:Nc
    rawDataTmp{i}=double(rawData((i-1)*sx+1:i*sx,:));
end

%% prepare phase correction to correct for offcenter data

phaseOff = phaseOffsetCorrection_BruKitchen(bruk);
%% sort data according to RARE factore/ Nb of echo / Nb of channel

im2=zeros(sx,sy,sz,NEchos,Nc,NR);
imOut=zeros(sx,sy,sz,NEchos,NR);
kdata=zeros(sx,sy,sz,NEchos,NR);

for k=1:NR
    for i=1:Nc
        temp=reshape(rawDataTmp{i}(1:Genpoint*NEchos*sx*NR),sx,RareFactor,NEchos,[],NR);
        for j=1:NEchos
            temp2=reshape(temp(:,:,j,:,k),sx,[]);
            rawData=zeros(sx,sy,sz);
            for yy=1:Genpoint
                rawData(:,GradPhaseVector(yy),GradSliceVector(yy))=temp2(:,yy);
            end
            
            rawData=reshape(rawData,sx,sy,sz).*phaseOff;
            kdata(:,:,:,j,i,k)=squeeze(double(rawData));
            
            im2(:,:,:,j,i,k)=fftshift(ifft(ifft(ifft(ifftshift(double(rawData)),[],1),[],2),[],3));
        end
    end
    imOut(:,:,:,:,k)=sqrt(sum(abs(im2(:,:,:,:,:,k)).^2,5));
end

%% Create a dataset for bart
% dim 1/2/3/4 =RO, PE, SE, Ch
% dim 6 contrast
% dim 8 Repetition ; 
kdata_bart = permute(kdata, [1 2 3 5 7 4 8 6 ]); 
%kdata_bart = permute(kdata, [1 2 3 5 7 4 11 8 9 10 6]); % Put NR data along 11 dimension

if strcmp(s_in.OPTION_RECO,'K_AVERAGE')
    % kdata_bart = mean(kdata_bart,8);
    tmp = sum(kdata_bart~=0,8);
    tmp(tmp == 0) = 1;
    kdata_bart = sum(kdata_bart,8) ./ tmp;
end

calibsz = (bruk.method.CenterMaskSize); % mask central fait 20*10 à cause de la résolution/FOV (128*128

if s_in.GPU
    warning('GPU is used : if reconstruction failder try with GPU = 0');
    
    sens_bart = bart(sprintf('ecalib -g -S -c0.1 -r %d -m %d',calibsz,s_in.SENS_MAPS), kdata_bart(:,:,:,:,:,2,1)); % second TI and first rep is used
    imCS = bart(sprintf('pics -d5 -g -U -i%d -e -R W:7:0:%d -R T:7:0:%d',s_in.ITER,s_in.W_LAMBDA,s_in.TV_LAMBDA),  kdata_bart , sens_bart);
else
    sens_bart = bart(sprintf('ecalib -S -c0.1 -r %d -m %d',calibsz,s_in.SENS_MAPS), kdata_bart(:,:,:,:,:,2,1)); % second TI and first rep is used
    imCS = bart(sprintf('pics -d5 -U -i%d -e -R W:7:0:%d -R T:7:0:%d',s_in.ITER,s_in.W_LAMBDA,s_in.TV_LAMBDA),  kdata_bart , sens_bart);
end

if strcmp(s_in.OPTION_RECO,'IM_AVERAGE')
    imCS=squeeze(mean(imCS(:,:,:,:,1,:,:,:),8)); % use only first set of map
end

%% Squeeze data to obtain : RO, PE,SE,TI,NR
imCS = squeeze(imCS);
%% MP2RAGE combinaison
    
struct_MP2RAGE.calibSize = calibsz;
struct_MP2RAGE.ETL = RareFactor;
struct_MP2RAGE.TI1 = bruk.method.EffectiveTE(1);
struct_MP2RAGE.TI2 = bruk.method.EffectiveTE(2);
struct_MP2RAGE.alpha1 = bruk.method.ExcPulse1(3);
struct_MP2RAGE.alpha2 = bruk.method.ExcPulse2(3);
struct_MP2RAGE.MP2RAGE_TR = bruk.method.RecoveryTime;
struct_MP2RAGE.TR = bruk.method.PVM_RepetitionTime;


[s_out] = MP2RAGE_lookupTable(imCS,struct_MP2RAGE); %% comb then MP2RAGE
disp('-----------RECO WITH correction background---------------');

bg_mult = 1;
multiFactor=bg_mult*mean(mean(mean(abs(imCS(1:end,end-10:end,end-10:end,2)))));
MP2RAGE_mask=real((conj(imCS(:,:,:,1,:)).*imCS(:,:,:,2,:)-multiFactor)./(abs(imCS(:,:,:,1,:)).^2+abs(imCS(:,:,:,2,:)).^2+2*multiFactor));

%% struct out

s_out.MP2RAGE_mask = squeeze(MP2RAGE_mask);
s_out.imTI = imCS;
s_out.param_in = s_in;
toc
end
