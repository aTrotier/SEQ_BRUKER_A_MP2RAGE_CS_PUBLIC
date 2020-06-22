%% function for simulation of SNR

%% Set parameters for simulation : stored in structure s

% choose 2 T1 tissues to optimize MP2RAGE contrast between them :
T1_WM = 1500
T1_GM = 900


% param MP2RAGE
ETL_list = [1:5:200];

% our parameters
TI1 = 400:100:1500;
TI2 = 1500:100:3000;
MP2RAGE_TR = 4000:100:8000;
alpha1 = 1:12;
alpha2 = 1:12;
TR = 7;


%% Calcul MP2RAGE
s.TR = TR;
% boucle TI1/TI2/MP2RAGE_TR/ETL/ALPHA1/ALPHA2
CNR = zeros(length(TI1),length(TI2),length(MP2RAGE_TR),length(ETL_list),length(alpha1),length(alpha2));
CNR2 = zeros(length(TI1),length(TI2),length(MP2RAGE_TR),length(ETL_list),length(alpha1),length(alpha2));
for idx_TI1 = 1:length(TI1)
    s.TI1 = TI1(idx_TI1);
    
    for idx_TI2 = 1:length(TI2)
        s.TI2 = TI2(idx_TI2);
        
        for idx_MP2RAGE_TR = 1:length(MP2RAGE_TR)
            s.MP2RAGE_TR = MP2RAGE_TR(idx_MP2RAGE_TR);
            
            for idx_ETL=1:length(ETL_list)
                s.ETL = ETL_list(idx_ETL);
                
                for idx_alpha1 = 1:length(alpha1)
                    s.alpha1 = alpha1(idx_alpha1);
                    
                    for idx_alpha2 = 1:length(alpha2)
                        s.alpha2 = alpha2(idx_alpha2);
                        
                        [S_MP2RAGE_WM, N_MP2RAGE_WM] = signal_MP2RAGE(s,T1_WM);
                        [S_MP2RAGE_GM, N_MP2RAGE_GM] = signal_MP2RAGE(s,T1_GM);
                        
                        
                        CNR(idx_TI1,idx_TI2,idx_MP2RAGE_TR,idx_ETL,idx_alpha1,idx_alpha2) = (S_MP2RAGE_WM - S_MP2RAGE_GM)/sqrt(N_MP2RAGE_WM^2 + N_MP2RAGE_GM^2)*1/sqrt(s.MP2RAGE_TR/1000);
                        CNR2(idx_TI1,idx_TI2,idx_MP2RAGE_TR,idx_ETL,idx_alpha1,idx_alpha2) = (S_MP2RAGE_WM - S_MP2RAGE_GM)/sqrt(N_MP2RAGE_WM^2 + N_MP2RAGE_GM^2)*sqrt(s.ETL)*1/sqrt(s.MP2RAGE_TR/1000);
                    end
                end
            end
        end
    end
end

[I,C] = max(abs(CNR2),[],'all','linear')

sz = [length(TI1),length(TI2),length(MP2RAGE_TR),length(ETL_list),length(alpha1),length(alpha2)];
[idx1,idx2,idx3,idx4,idx5,idx6] = ind2sub(sz,C);
% TI1/TI2/MP2RAGE_TR/ETL/ALPHA1/ALPHA2
disp(['CNR2 -> TI1=' num2str(TI1(idx1)) ' TI2=' num2str(TI2(idx2)) ' MP2RAGE_TR=' num2str(MP2RAGE_TR(idx3)) ' ETL=' num2str(ETL_list(idx4)) ' alpha1=' num2str(alpha1(idx5)) ' alpha2=' num2str(alpha2(idx6)) ]);

%
[I,C] = max(abs(CNR),[],'all','linear')

sz = [length(TI1),length(TI2),length(MP2RAGE_TR),length(ETL_list),length(alpha1),length(alpha2)];
[idx1,idx2,idx3,idx4,idx5,idx6] = ind2sub(sz,C);
% TI1/TI2/MP2RAGE_TR/ETL/ALPHA1/ALPHA2
disp(['CNR -> TI1=' num2str(TI1(idx1)) ' TI2=' num2str(TI2(idx2)) ' MP2RAGE_TR=' num2str(MP2RAGE_TR(idx3)) ' ETL=' num2str(ETL_list(idx4)) ' alpha1=' num2str(alpha1(idx5)) ' alpha2=' num2str(alpha2(idx6)) ]);

%% Function
function [S_MP2RAGE, N_MP2RAGE] = signal_MP2RAGE(struct_MP2RAGE,T1)
n = struct_MP2RAGE.ETL;
TI1 = struct_MP2RAGE.TI1;
TI2 = struct_MP2RAGE.TI2;
alpha1 = struct_MP2RAGE.alpha1;
alpha2 = struct_MP2RAGE.alpha2;
MP2RAGETR = struct_MP2RAGE.MP2RAGE_TR;
TR = struct_MP2RAGE.TR;

% check parameter
% disp('-------------------------------------------------');
% disp(['alpha1 :', num2str(alpha1), ' alpha2 :', num2str(alpha2)]);
% disp(['TI1 :', num2str(TI1), ' TI2 :', num2str(TI2)]);
% disp(['TR :', num2str(TR), ' MP2RAGETR :', num2str(MP2RAGETR), ' size of echo trains :', num2str(n)]);
% disp('-------------------------------------------------');

E1=exp(-TR./T1);
EA=exp(-(TI1-(n./2-1).*TR)./T1);
EB=exp(-(TI2-TI1-n.*TR)./T1);
EC=exp(-(MP2RAGETR-(TI2+(n./2).*TR))./T1);

eff=1;

% compute mzss =[(B).*(cosd(alpha2).*E1).^n+A].*EC+(1-EC);

B = ((1-EA).*(cosd(alpha1).*E1).^n+(1-E1).*(1-(cosd(alpha1).*E1).^n)./(1-cosd(alpha1).*E1)).*EB+(1-EB);
A = (1-E1).*((1-(cosd(alpha2).*E1).^n)./(1-cosd(alpha2).*E1));

mzss_num=((B).*(cosd(alpha2).*E1).^n+A).*EC+(1-EC);
mzss_denom=(1+eff.*(cosd(alpha1).*cosd(alpha2)).^n .* exp(-MP2RAGETR./T1));

mzss=mzss_num./mzss_denom;

% compute GRE1= sind(alpha1).*(A.*(cosd(alpha1).*E1).^(n./2-1)+B)

A=-eff.*mzss.*EA+(1-EA);
B=(1-E1).*(1-(cosd(alpha1).*E1).^(n./2-1))./(1-cosd(alpha1).*E1);

GRE1=sind(alpha1).*(A.*(cosd(alpha1).*E1).^(n./2-1)+B);

% compute GRE2= sind(alpha2).*(A-B)

A=(mzss-(1-EC))./(EC.*(cosd(alpha2).*E1).^(n./2));
B=(1-E1).*((cosd(alpha2).*E1).^(-n./2)-1)./(1-cosd(alpha2).*E1);

GRE2=sind(alpha2).*(A-B);

S_MP2RAGE=GRE2.*GRE1./(GRE1.*GRE1+GRE2.*GRE2);


% noise of the sequence :

N_MP2RAGE = sqrt((GRE1^2-GRE2^2)^2/(GRE1^2+GRE2^2)^3);
end
