function [s_out] = reco_multi_param_W(s_in)
s_in.kdata_fully;
s_in.kdata_under;

s_in.Wlist;
s_in.ITERlist;

s_in.GPU;
s_in.calib;
s_in.SENS_MAPS;

%% sensitivity maps :

if s_in.GPU
    sens_bart = bart(sprintf('ecalib -g -d3  -r %d -m %d',s_in.calib,s_in.SENS_MAPS), s_in.kdata_under(:,:,:,:,:,2,1)); % second TI and first rep is used
else
    sens_bart = bart(sprintf('ecalib -d3  -r %d -m %d',s_in.calib,s_in.SENS_MAPS), s_in.kdata_under(:,:,:,:,:,2,1)); % second TI and first rep is used
end

%% Fully reco : only use the first map

im_fft = bart('pics -d4 -i1 -e',s_in.kdata_fully,sens_bart);
im_fft = abs(squeeze(im_fft(:,:,:,:,1,2))); % only first set of maps used and 2nd TI
%% Zero filling

im_zf = bart('pics -d4 -i1 -e',s_in.kdata_under,sens_bart);
im_zf = abs(squeeze(im_zf(:,:,:,:,1,2))); % only first set of maps used and 2nd TI


%% LOOP

tmp1 = 10000;
for w = 1:length(s_in.Wlist)
    for it = 1:length(s_in.ITERlist)
        fprintf(' W = %d - iter = %d',s_in.Wlist(w),s_in.ITERlist(it));
        
        if s_in.GPU
            im_CS = bart(sprintf('pics -d3 -g -U -i%d -e -R W:7:0:%d',s_in.ITERlist(it),s_in.Wlist(w)),  s_in.kdata_under , sens_bart);
            %im_CS2 = bart(sprintf('pics -d3 -g -U -i%d -e -R W:7:0:%d -R W:32:0:%d',s_in.ITERlist(it),s_in.Wlist(w),s_in.Wlist2(tv)),  s_in.kdata_under , sens_bart);
        else
            im_CS = bart(sprintf('pics -d3 -U -i%d -e -R W:7:0:%d',s_in.ITERlist(it),s_in.Wlist(w)),  s_in.kdata_under , sens_bart);
            %im_CS2 = bart(sprintf('pics -d3 -U -i%d -e -R W:7:0:%d -R W:32:0:%d',s_in.ITERlist(it),s_in.Wlist(w),s_in.Wlist2(tv)),  s_in.kdata_under , sens_bart);
        end
        
        im_CS = abs(squeeze(im_CS(:,:,:,:,1,2))); % only first set of maps used and 2nd TI
        %% Comparaison
        diff_CS = (im_CS-im_fft);
        s_out.RMSD(w,it) = rms(diff_CS(:));
        
        if(s_out.RMSD(w,it) < tmp1)
            s_out.im_CS = im_CS;
            tmp1 = s_out.RMSD(w,it);
        end
        
        
    end
end



s_out.im_fft = im_fft;
s_out.im_zf = im_zf;

s_out.diff_CS = (s_out.im_CS-s_out.im_fft);

s_out.diff_zf = (s_out.im_zf-s_out.im_fft);


s_out.RMSD_zf = rms(s_out.diff_zf(:));
fprintf('RMSD of zeroFilling = %d', s_out.RMSD_zf)
end