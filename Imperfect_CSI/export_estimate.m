
          
% par.detector = {'Gauss-Seidel1','Newton_richardson', 'Newton_chebyshev','Newton_GS',...
%     'Conjugate-Gradient', 'richardson', 'MMSE', 'chebyshev',  ...
%     'LAS-bMMSE', 'LAS-Newton_richardson', 'LAS-Newton_chebyshev', 'LAS-Newton_GS',...
%     'LAS-GS', 'LAS-richardson','LAS-chebyshev'};          


fullname = 'results/estimate_BER_No_correlation32x128_64QAM_RunID_0.mat';
%  {'Gauss-Seidel1','Newton_richardson', 'Newton_chebyshev','Newton_GS',...
%                'Conjugate-Gradient', 'richardson', 'bMMSE', 'chebyshev','steepest_jacobi','MMSEapprox'}
tmp = load(fullname);
SNR                 = tmp.Results.SNRdB_list';
GS                  = tmp.Results.BER(1,:)';
NS_RI               = tmp.Results.BER(2,:)';
NS_Cheby            = tmp.Results.BER(3,:)';
NS_GS               = tmp.Results.BER(4,:)';
CG                  = tmp.Results.BER(5,:)';
richardson          = tmp.Results.BER(6,:)';
bMMSE               = tmp.Results.BER(7,:)';
Cheby               = tmp.Results.BER(8,:)';      
LAS_MMSE            = tmp.Results.BER(9,:)';
LAS_NS_RI           = tmp.Results.BER(10,:)';
LAS_NS_Cheby        = tmp.Results.BER(11,:)';
LAS_NS_GS           = tmp.Results.BER(12,:)';
LAS_GS              = tmp.Results.BER(13,:)';
LAS_RI              = tmp.Results.BER(14,:)';
LAS_Cheby           = tmp.Results.BER(15,:)';
 
 
 
FileName = 'Compiled_Data_Origin/estimate_64QAM_32x128.mat';
save(FileName,'SNR','GS', 'NS_RI',...
              'NS_Cheby', 'NS_GS',... 
              'CG', 'richardson','bMMSE','Cheby','LAS_MMSE',...
              'LAS_NS_RI', 'LAS_NS_Cheby','LAS_NS_GS', 'LAS_GS','LAS_RI', 'LAS_Cheby');
          
          

              
fullname = 'results/estimate_BER_moderate_025correlation32x128_64QAM_RunID_0.mat';
 
tmp = load(fullname);
SNR                 = tmp.Results.SNRdB_list';
GS                  = tmp.Results.BER(1,:)';
NS_RI               = tmp.Results.BER(2,:)';
NS_Cheby            = tmp.Results.BER(3,:)';
NS_GS               = tmp.Results.BER(4,:)';
CG                  = tmp.Results.BER(5,:)';
richardson          = tmp.Results.BER(6,:)';
bMMSE               = tmp.Results.BER(7,:)';
Cheby               = tmp.Results.BER(8,:)';      
LAS_MMSE            = tmp.Results.BER(9,:)';
LAS_NS_RI           = tmp.Results.BER(10,:)';
LAS_NS_Cheby        = tmp.Results.BER(11,:)';
LAS_NS_GS           = tmp.Results.BER(12,:)';
LAS_GS              = tmp.Results.BER(13,:)';
LAS_RI              = tmp.Results.BER(14,:)';
LAS_Cheby           = tmp.Results.BER(15,:)';
 
 
 
FileName = 'Compiled_Data_Origin/estimate_modarate025_64QAM_32x128.mat';
save(FileName,'SNR','GS', 'NS_RI',...
              'NS_Cheby', 'NS_GS',... 
              'CG', 'richardson','bMMSE','Cheby','LAS_MMSE',...
              'LAS_NS_RI', 'LAS_NS_Cheby','LAS_NS_GS', 'LAS_GS','LAS_RI', 'LAS_Cheby');
