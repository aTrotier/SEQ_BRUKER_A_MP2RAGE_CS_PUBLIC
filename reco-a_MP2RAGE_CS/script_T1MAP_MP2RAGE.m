% Script illustrating how to reconstruct a Bruker dataset acquired with the
% sequence a_MP2RAGE_CS
%
% TO DO : add nifti export


%% OPTIONAL CHANGE OMP_NUM_THREADS in order to accelerate BART reconstruction
% num_PhysCPUs = feature('numcores')
% setenv('OMP_NUM_THREADS', num2str(num_PhysCPUs));
%% We use an object to create the parameter for the reconstruction 

param_in = OBJ_MP2RAGE_RECO; %Without input it will ask for the brucker dataset folder (where fid is located)
% we can also directly give the path if wanted with :
% param_in = OBJ_MP2RAGE_RECO(path);

%% Edit the parameter for the reconstruction

% param_in = 
% 
%   OBJ_MP2RAGE_RECO with properties:
% 
%     BRUKER_PATH: '/home/atrotier/GITHUB/SEQ_BRUKER_a_MP2RAGE_CS/data_test/'
%             GPU: 0
%     OPTION_RECO: "K_AVERAGE"
%            ITER: 60
%       SENS_MAPS: 1
%        W_LAMBDA: 0.0100
%       TV_LAMBDA: 1.0000e-03

param_in.ITER = 30; % iteration for 
param_in.OPTION_RECO = "NR_RECO"; % reconstruction of NR is done independently 
% other possibilities are : "K_AVERAGE" and "IM_AVERAGE"

%% Run the reconstruction

s_out = a_MP2RAGE_CS_bart(param_in);

figure; imagesc(s_out.T1map(:,:,end/2)); colormap('gray'); colorbar;

