%
% imOut=phaseOffsetCorrection()
%
% Author:   Aur�lien TROTIER  (a.trotier@gmail.com)
% Date:     2016-06-22
% Partner: none
% Institute: CRMSB (Bordeaux)
%
% Function description:

%
% Input:
%       bruker structure : generate using BruKitchen
%
% Output:
%
%   phaseOff = exp(-i*2pi*phase)  à multiplier avec le kspace
%
%
% Algorithm & bibliography:
%
% See also :
%
% To do :
%   *fonctionne pour la 3D à modifier pour la 2D

function phaseOff = phaseOffsetCorrection_BruKitchen(bruk)

sx=bruk.acqp.ACQ_size(1)/2;
sy=bruk.method.PVM_Matrix(2);
PVM_SpatDimEnum = bruk.method.PVM_SpatDimEnum{1};


if strcmp(PVM_SpatDimEnum,'3D')
    sz=bruk.method.PVM_Matrix(3);
    phaseOffy=zeros(sx,sy,sz);
    phaseOffz=zeros(sx,sy,sz);
    
    for k=1:sz
        phaseOffz(:,:,k)=ones(sx,sy)*(k-sz/2+1)/sz*sz*bruk.method.PVM_SPackArrSliceOffset/bruk.method.PVM_Fov(3);
    end
else
    sz=1;
    phaseOffy=zeros(sx,sy,sz);
    phaseOffz=zeros(sx,sy,sz);
end




for k=1:sy
    phaseOffy(:,k,:)=ones(sx,sz)*(k-sy/2+1)/sy*sy*bruk.method.PVM_SPackArrPhase1Offset/bruk.method.PVM_Fov(2);
end

phaseOff=phaseOffy+phaseOffz;


phaseOff=exp(1i*2*pi*phaseOff);

end