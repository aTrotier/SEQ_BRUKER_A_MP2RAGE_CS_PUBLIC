classdef OBJ_MP2RAGE_RECO
    %OBJ_MP2RAGE_RECO class used as input for the function a_MP2RAGE_CS_bart
    %   
    
    properties
        BRUKER_PATH char = "";
        
        GPU(1,1) {mustBeMember(GPU,[0 1])} = 0;
        OPTION_RECO(1,1) string {mustBeMember(OPTION_RECO, ["K_AVERAGE","IM_AVERAGE","NR_RECO"])} = "K_AVERAGE";
        
        ITER(1,1) {mustBeInteger, mustBePositive} = 60;
        SENS_MAPS(1,1) {mustBeInteger, mustBePositive} = 1;
        W_LAMBDA(1,1) {mustBeNumeric, mustBeNonnegative} = 0.01;
        TV_LAMBDA(1,1) {mustBeNumeric, mustBeNonnegative} = 0.001;
    end
    methods
        %% constructor
        function obj = OBJ_MP2RAGE_RECO(path)
            if nargin < 1
                path = uigetdir('*.*','Open rawdata bruker directory');
                cd(path);
            end
            
            if ischar(path)
                obj.BRUKER_PATH = path;
            else
                error("input is not a path to a bruker dataset folder");
            end
        end
        %% methods

    end
end

