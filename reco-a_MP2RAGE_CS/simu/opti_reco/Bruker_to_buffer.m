
function kdata_bart = Bruker_to_buffer(bruk)
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

% read fid file in the same directory
rawData=bruk.fid;
% Remove 0 from data
sizeR2=128*ceil(sx*Nc/128);
rawData=reshape(rawData,sizeR2,[]);
rawData=rawData(1:Nc*sx,:);

rawDataTmp=cell(Nc,1);
for i=1:Nc
    rawDataTmp{i}=double(rawData((i-1)*sx+1:i*sx,:));
end

% prepare phase correction to correct for offcenter data

phaseOff = phaseOffsetCorrection_BruKitchen(bruk);

% sort data according to RARE factore/ Nb of echo / Nb of channel

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

% Create a dataset for bart
% dim 1/2/3/4 =RO, PE, SE, Ch
% dim 6 contrast
% dim 8 Repetition ;
kdata_bart = permute(kdata, [1 2 3 5 7 4 8 6 ]);

end
