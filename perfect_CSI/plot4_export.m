
%% only MMSE section
clear
fullname = 'results/Plot4_BER_32x128_256QAM_RunID_0.mat';
%  {'Gauss-Seidel1','Newton_richardson', 'Newton_chebyshev','Newton_GS',...
%                'Conjugate-Gradient', 'richardson', 'bMMSE', 'chebyshev','steepest_jacobi','MMSEapprox'}
tmp = load(fullname);
SNR                      = tmp.Results.SNRdB_list';
GS256QAM                  = tmp.Results.BER(1,:)';
Newton_richardson256QAM   = tmp.Results.BER(2,:)';
Newton_chebyshev256QAM    = tmp.Results.BER(3,:)';
Newton_GS256QAM           = tmp.Results.BER(4,:)';
CG256QAM                  = tmp.Results.BER(5,:)';
richardson256QAM          = tmp.Results.BER(6,:)';
bMMSE                    = tmp.Results.BER(7,:)';
chebyshev256QAM           = tmp.Results.BER(8,:)';
steepest_jacobi256QAM     = tmp.Results.BER(9,:)';
MMSEapprox256QAM          = tmp.Results.BER(10,:)';
 
 
 
 
FileName = 'Compiled_Data_Origin/plot4_256QAM_32x128.mat';
save(FileName,'SNR','GS256QAM', 'Newton_richardson256QAM',...
              'Newton_chebyshev256QAM', 'Newton_GS256QAM',... 
              'CG256QAM', 'richardson256QAM','bMMSE','chebyshev256QAM',...
              'steepest_jacobi256QAM', 'MMSEapprox256QAM');
          
          


fullname = 'results/Plot4_BER_32x128_64QAM_RunID_0.mat';
%  {'Gauss-Seidel1','Newton_richardson', 'Newton_chebyshev','Newton_GS',...
%                'Conjugate-Gradient', 'richardson', 'bMMSE', 'chebyshev','steepest_jacobi','MMSEapprox'}
tmp = load(fullname);
SNR                      = tmp.Results.SNRdB_list';
GS64QAM                  = tmp.Results.BER(1,:)';
Newton_richardson64QAM   = tmp.Results.BER(2,:)';
Newton_chebyshev64QAM    = tmp.Results.BER(3,:)';
Newton_GS64QAM           = tmp.Results.BER(4,:)';
CG64QAM                  = tmp.Results.BER(5,:)';
richardson64QAM          = tmp.Results.BER(6,:)';
bMMSE                    = tmp.Results.BER(7,:)';
chebyshev64QAM           = tmp.Results.BER(8,:)';
steepest_jacobi64QAM     = tmp.Results.BER(9,:)';
MMSEapprox64QAM          = tmp.Results.BER(10,:)';
 
 
 
 
FileName = 'Compiled_Data_Origin/plot4_64QAM_32x128.mat';
save(FileName,'SNR','GS64QAM', 'Newton_richardson64QAM',...
              'Newton_chebyshev64QAM', 'Newton_GS64QAM',... 
              'CG64QAM', 'richardson64QAM','bMMSE','chebyshev64QAM',...
              'steepest_jacobi64QAM', 'MMSEapprox64QAM');
          
          

              
fullname = 'results/Plot4_BER_32x128_16QAM_RunID_0.mat';
 
tmp = load(fullname);
SNR                      = tmp.Results.SNRdB_list';
GS16QAM                  = tmp.Results.BER(1,:)';
Newton_richardson16QAM   = tmp.Results.BER(2,:)';
Newton_chebyshev16QAM    = tmp.Results.BER(3,:)';
Newton_GS16QAM           = tmp.Results.BER(4,:)';
CG16QAM                  = tmp.Results.BER(5,:)';
richardson16QAM          = tmp.Results.BER(6,:)';
bMMSE                    = tmp.Results.BER(7,:)';
chebyshev16QAM           = tmp.Results.BER(8,:)';
steepest_jacobi16QAM     = tmp.Results.BER(9,:)';
MMSEapprox16QAM          = tmp.Results.BER(10,:)';
 
 
 
 
FileName = 'Compiled_Data_Origin/plot4_16QAM_32x128.mat';
save(FileName,'SNR','GS16QAM', 'Newton_richardson16QAM',...
              'Newton_chebyshev16QAM', 'Newton_GS16QAM',... 
              'CG16QAM', 'richardson16QAM','bMMSE','chebyshev16QAM',...
              'steepest_jacobi16QAM', 'MMSEapprox16QAM');
          
          
%% LAS section 
% 64QAM

fullname1 = 'results/Plot4_BER_LAS_full_32x128_64QAM_RunID_0.mat';
fullname2 = 'results/Plot4_BER_LAS_partial_w_2032x128_64QAM_RunID_0.mat';
 
tmp = load(fullname1);
SNR                            = tmp.Results.SNRdB_list';
bMMSE64QAM                     = tmp.Results.BER(1,:)';
LAS_bMMSE64QAM_f               = tmp.Results.BER(2,:)';
LAS_Newton_richardson64QAM_f   = tmp.Results.BER(3,:)';
LAS_Newton_chebyshev64QAM_f    = tmp.Results.BER(4,:)';
LAS_Newton_GS64QAM_f           = tmp.Results.BER(5,:)';
LAS_GS64QAM_f                  = tmp.Results.BER(6,:)';
LAS_richardson64QAM_f          = tmp.Results.BER(7,:)';
LAS_chebyshev64QAM_f           = tmp.Results.BER(8,:)';


 
 
 
FileName = 'Compiled_Data_Origin/plot4_LAS_64QAM_32x128.mat';
save(FileName,'SNR','bMMSE64QAM','LAS_bMMSE64QAM_f', ...
    'LAS_Newton_richardson64QAM_f','LAS_Newton_chebyshev64QAM_f',...
    'LAS_Newton_GS64QAM_f', 'LAS_GS64QAM_f',  ...
    'LAS_richardson64QAM_f','LAS_chebyshev64QAM_f');
                    
tmp2 = load(fullname2);
bMMSE64QAM                     = tmp.Results.BER(1,:)';
LAS_bMMSE64QAM_p               = tmp2.Results.BER(2,:)';
LAS_Newton_richardson64QAM_p   = tmp2.Results.BER(3,:)';
LAS_Newton_chebyshev64QAM_p    = tmp2.Results.BER(4,:)';
LAS_Newton_GS64QAM_p           = tmp2.Results.BER(5,:)';
LAS_GS64QAM_p                  = tmp2.Results.BER(6,:)';
LAS_richardson64QAM_p          = tmp2.Results.BER(7,:)';
LAS_chebyshev64QAM_p           = tmp2.Results.BER(8,:)'; 


FileName = 'Compiled_Data_Origin/plot4_partial_LAS_64QAM_32x128.mat';
save(FileName,'SNR','bMMSE64QAM','LAS_bMMSE64QAM_p',...
    'LAS_Newton_richardson64QAM_p','LAS_Newton_chebyshev64QAM_p',...
    'LAS_Newton_GS64QAM_p', 'LAS_GS64QAM_p', ...
    'LAS_richardson64QAM_p','LAS_chebyshev64QAM_p');

% 16QAM

% fullname1 = 'results/Plot4_BER_LAS_full_32x128_16QAM_RunID_0.mat';
% fullname2 = 'results/Plot4_BER_LAS_partial_w_8_32x128_16QAM_RunID_0.mat';
%  
% tmp = load(fullname1);
% SNR                            = tmp.Results.SNRdB_list';
% bMMSE16QAM                     = tmp.Results.BER(1,:)';
% LAS_bMMSE16QAM_f               = tmp.Results.BER(2,:)';
% LAS_Newton_richardson16QAM_f   = tmp.Results.BER(3,:)';
% LAS_Newton_chebyshev16QAM_f    = tmp.Results.BER(4,:)';
% LAS_Newton_GS16QAM_f           = tmp.Results.BER(5,:)';
% LAS_GS16QAM_f                  = tmp.Results.BER(6,:)';
% LAS_richardson16QAM_f          = tmp.Results.BER(7,:)';
% LAS_chebyshev16QAM_f           = tmp.Results.BER(8,:)';
% 
% tmp2 = load(fullname2);
% LAS_bMMSE16QAM_p               = tmp2.Results.BER(2,:)';
% LAS_Newton_richardson16QAM_p   = tmp2.Results.BER(3,:)';
% LAS_Newton_chebyshev16QAM_p    = tmp2.Results.BER(4,:)';
% LAS_Newton_GS16QAM_p           = tmp2.Results.BER(5,:)';
% LAS_GS16QAM_p                  = tmp2.Results.BER(6,:)';
% LAS_richardson16QAM_p          = tmp2.Results.BER(7,:)';
% LAS_chebyshev16QAM_p           = tmp2.Results.BER(8,:)'; 
%  
%  
%  
% FileName = 'Compiled_Data_Origin/plot4_LAS_16QAM_32x128.mat';
% save(FileName,'SNR','bMMSE16QAM','LAS_bMMSE16QAM_f', 'LAS_bMMSE16QAM_p',...
%     'LAS_Newton_richardson16QAM_f', 'LAS_Newton_richardson16QAM_p',...
%     'LAS_Newton_chebyshev16QAM_f', 'LAS_Newton_chebyshev16QAM_p',...
%     'LAS_Newton_GS16QAM_f', 'LAS_Newton_GS16QAM_p', ...
%     'LAS_GS16QAM_f', 'LAS_GS16QAM_p', ...
%     'LAS_richardson16QAM_f', 'LAS_richardson16QAM_p', ...
%     'LAS_chebyshev16QAM_f', 'LAS_chebyshev16QAM_p');
% 
%  
