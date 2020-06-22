function [s_out] = reco_multi_param(s_in)
 s_in.kdata_fully;
 s_in.kdata_under;
 
 s_in.Wlist;
 s_in.TVlist;
 s_in.ITERlist;

 s_in.GPU;
 s_in.calib;
 s_in.SENS_MAPS;
 
 bg_mult=s_in.bg_mult;
 
 %% sensitivity maps :

if s_in.GPU
    sens_bart = bart(sprintf('ecalib -d3 -c0 -g -r %d -m %d',s_in.calib,s_in.SENS_MAPS), s_in.kdata_under(:,:,:,:,:,2,1)); % second TI and first rep is used
else
    sens_bart = bart(sprintf('ecalib -d3 -c0 -r %d -m %d',s_in.calib,s_in.SENS_MAPS), s_in.kdata_under(:,:,:,:,:,2,1)); % second TI and first rep is used
end

 %% Fully reco : only use the first map

im_fft = bart('pics -d4 -i1 -e',s_in.kdata_fully,sens_bart);
im_fft = abs(squeeze(im_fft(:,:,:,:,1,:))); % only first set of maps used but need the 2 TI for MP2RAGE


multiFactor=bg_mult*mean(mean(mean(abs(im_fft(1:end,end-10:end,end-10:end,2)))));
MP2RAGE_fft=real((conj(im_fft(:,:,:,1,:)).*im_fft(:,:,:,2,:)-multiFactor)./(abs(im_fft(:,:,:,1,:)).^2+abs(im_fft(:,:,:,2,:)).^2+2*multiFactor));

%% Zero filling

im_zf = bart('pics -d4 -i1 -e',s_in.kdata_under,sens_bart);
im_zf = abs(squeeze(im_zf(:,:,:,:,1,:))); % only first set of maps used and 2nd TI


%multiFactor=bg_mult*mean(mean(mean(abs(im_zf(1:end,end-10:end,end-10:end,2)))));
MP2RAGE_zf=real((conj(im_zf(:,:,:,1,:)).*im_zf(:,:,:,2,:)-multiFactor)./(abs(im_zf(:,:,:,1,:)).^2+abs(im_zf(:,:,:,2,:)).^2+2*multiFactor));

%% LOOP 
tmp2 = 10000;
tmp1 = 10000;
tmp3 = 0;
tmp4 = 0;
for tv = 1:length(s_in.TVlist)
    for w = 1:length(s_in.Wlist)
        for it = 1:length(s_in.ITERlist)
            fprintf('TV = %d - W = %d - iter = %d',s_in.TVlist(tv),s_in.Wlist(w),s_in.ITERlist(it));
            
            if s_in.GPU
                im_CS = bart(sprintf('pics -d3 -g -U -i%d -e -R W:7:0:%d -R T:7:0:%d',s_in.ITERlist(it),s_in.Wlist(w),s_in.TVlist(tv)),  s_in.kdata_under , sens_bart);
                im_CS2 = bart(sprintf('pics -d3 -g -U -i%d -e -R W:7:0:%d -R T:32:0:%d',s_in.ITERlist(it),s_in.Wlist(w),s_in.TVlist(tv)),  s_in.kdata_under , sens_bart);
            else
                im_CS = bart(sprintf('pics -d3 -U -i%d -e -R W:7:0:%d -R T:7:0:%d',s_in.ITERlist(it),s_in.Wlist(w),s_in.TVlist(tv)),  s_in.kdata_under , sens_bart);
                im_CS2 = bart(sprintf('pics -d3 -U -i%d -e -R W:7:0:%d -R T:32:0:%d',s_in.ITERlist(it),s_in.Wlist(w),s_in.TVlist(tv)),  s_in.kdata_under , sens_bart);
            end
            
            im_CS = abs(squeeze(im_CS(:,:,:,:,1,:))); % only first set of maps used and 2nd TI
            im_CS2 = abs(squeeze(im_CS2(:,:,:,:,1,:))); % only first set of maps used and 2nd TI
            
            %% create MP2RAGE image
            %multiFactor=bg_mult*mean(mean(mean(abs(im_CS(1:end,end-10:end,end-10:end,2)))));
            MP2RAGE_CS=real((conj(im_CS(:,:,:,1,:)).*im_CS(:,:,:,2,:)-multiFactor)./(abs(im_CS(:,:,:,1,:)).^2+abs(im_CS(:,:,:,2,:)).^2+2*multiFactor));
            
            %multiFactor=bg_mult*mean(mean(mean(abs(im_CS2(1:end,end-10:end,end-10:end,2)))));
            MP2RAGE_CS2=real((conj(im_CS2(:,:,:,1,:)).*im_CS2(:,:,:,2,:)-multiFactor)./(abs(im_CS2(:,:,:,1,:)).^2+abs(im_CS2(:,:,:,2,:)).^2+2*multiFactor));


            %% Comparaison         
            diff_CS = (MP2RAGE_CS-MP2RAGE_fft);            
            s_out.RMSD(tv,w,it) = rms(diff_CS(:));
            
            if(s_out.RMSD(tv,w,it) < tmp1)
                s_out.im_CS = im_CS;
                s_out.MP2RAGE_CS = MP2RAGE_CS;
                tmp1 = s_out.RMSD(tv,w,it);
            end
            
            s_out.SSIM(tv,w,it) = ssim(MP2RAGE_CS,MP2RAGE_fft);
            if(s_out.SSIM(tv,w,it) > tmp3)
                s_out.im_CS_SSIM = im_CS;
                s_out.MP2RAGE_CS_SSIM = MP2RAGE_CS;
                tmp3 = s_out.SSIM(tv,w,it);
            end
            
            diff_CS2 = (MP2RAGE_CS2-MP2RAGE_fft);            
            s_out.RMSD2(tv,w,it) = rms(diff_CS2(:));
            
            if(s_out.RMSD2(tv,w,it) < tmp2)
                s_out.im_CS2 = im_CS2;
                s_out.MP2RAGE_CS2 = MP2RAGE_CS2;
                tmp2 = s_out.RMSD2(tv,w,it);
            end
            
            s_out.SSIM2(tv,w,it) = ssim(MP2RAGE_CS2,MP2RAGE_fft);
            
            if(s_out.SSIM2(tv,w,it) > tmp4)
                s_out.im_CS_SSIM2 = im_CS2;
                s_out.MP2RAGE_CS_SSIM2 = MP2RAGE_CS2;
                tmp4 = s_out.SSIM2(tv,w,it);
            end
        end
    end
end

s_out.sens_bart = sens_barts;

s_out.im_fft = im_fft;
s_out.im_zf = im_zf;


s_out.MP2RAGE_fft = MP2RAGE_fft;
s_out.MP2RAGE_zf = MP2RAGE_zf;

s_out.SSIM_zf = ssim(MP2RAGE_zf,MP2RAGE_fft);
end