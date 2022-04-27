
%% only MMSE section
clear
par.detector = {'Gauss-Seidel1','Newton_richardson', 'Newton_chebyshev','Newton_GS',...
                 'Conjugate-Gradient', 'richardson', 'bMMSE', 'chebyshev', ...
                'LAS-bMMSE', 'LAS-Newton_richardson', 'LAS-Newton_chebyshev', 'LAS-Newton_GS',... 
                'LAS-GS', 'LAS-richardson','LAS-chebyshev'};             

            
fullname = 'results/Plot5_BER_32x128_64QAM_RunID_0.mat';

tmp = load(fullname);
rho                   = tmp.Results.rho';
GS                    = tmp.Results.BER(1,:)';
Newton_richardson     = tmp.Results.BER(2,:)';
Newton_chebyshev      = tmp.Results.BER(3,:)';
Newton_GS             = tmp.Results.BER(4,:)';
CG                    = tmp.Results.BER(5,:)';
richardson            = tmp.Results.BER(6,:)';
bMMSE                 = tmp.Results.BER(7,:)';
chebyshev             = tmp.Results.BER(8,:)';
LAS_bMMSE             = tmp.Results.BER(9,:)';
LAS_Newton_richardson = tmp.Results.BER(10,:)';
LAS_Newton_chebyshev  = tmp.Results.BER(11,:)';
LAS_Newton_GS         = tmp.Results.BER(12,:)';
LAS_GS                = tmp.Results.BER(13,:)';
LAS_richardson        = tmp.Results.BER(14,:)';
LAS_chebyshev         = tmp.Results.BER(15,:)';


FileName = 'Compiled_Data_Origin/plot5_correlation_64QAM_32x128.mat';
save(FileName,'rho','GS', 'Newton_richardson',...
              'Newton_chebyshev', 'Newton_GS',... 
              'CG', 'richardson','bMMSE', 'chebyshev', 'LAS_bMMSE', 'LAS_Newton_richardson', ...
              'LAS_Newton_chebyshev', 'LAS_Newton_GS', 'LAS_GS', 'LAS_richardson', 'LAS_chebyshev');

          
fullname = 'results/Plot5_BER_W16_32x128_64QAM_RunID_0.mat';

tmp = load(fullname);
rho                   = tmp.Results.rho';
GS                    = tmp.Results.BER(1,:)';
Newton_richardson     = tmp.Results.BER(2,:)';
Newton_chebyshev      = tmp.Results.BER(3,:)';
Newton_GS             = tmp.Results.BER(4,:)';
CG                    = tmp.Results.BER(5,:)';
richardson            = tmp.Results.BER(6,:)';
bMMSE                 = tmp.Results.BER(7,:)';
chebyshev             = tmp.Results.BER(8,:)';
LAS_bMMSE             = tmp.Results.BER(9,:)';
LAS_Newton_richardson = tmp.Results.BER(10,:)';
LAS_Newton_chebyshev  = tmp.Results.BER(11,:)';
LAS_Newton_GS         = tmp.Results.BER(12,:)';
LAS_GS                = tmp.Results.BER(13,:)';
LAS_richardson        = tmp.Results.BER(14,:)';
LAS_chebyshev         = tmp.Results.BER(15,:)';


FileName = 'Compiled_Data_Origin/plot5_correlationW16_64QAM_32x128.mat';
save(FileName,'rho','GS', 'Newton_richardson',...
              'Newton_chebyshev', 'Newton_GS',... 
              'CG', 'richardson','bMMSE', 'chebyshev', 'LAS_bMMSE', 'LAS_Newton_richardson', ...
              'LAS_Newton_chebyshev', 'LAS_Newton_GS', 'LAS_GS', 'LAS_richardson', 'LAS_chebyshev');          
    
fullname = 'results/Plot5_BER_W20_32x128_64QAM_RunID_0.mat';

tmp = load(fullname);
rho                   = tmp.Results.rho';
GS                    = tmp.Results.BER(1,:)';
Newton_richardson     = tmp.Results.BER(2,:)';
Newton_chebyshev      = tmp.Results.BER(3,:)';
Newton_GS             = tmp.Results.BER(4,:)';
CG                    = tmp.Results.BER(5,:)';
richardson            = tmp.Results.BER(6,:)';
bMMSE                 = tmp.Results.BER(7,:)';
chebyshev             = tmp.Results.BER(8,:)';
LAS_bMMSE             = tmp.Results.BER(9,:)';
LAS_Newton_richardson = tmp.Results.BER(10,:)';
LAS_Newton_chebyshev  = tmp.Results.BER(11,:)';
LAS_Newton_GS         = tmp.Results.BER(12,:)';
LAS_GS                = tmp.Results.BER(13,:)';
LAS_richardson        = tmp.Results.BER(14,:)';
LAS_chebyshev         = tmp.Results.BER(15,:)';


FileName = 'Compiled_Data_Origin/plot5_correlationW20_64QAM_32x128.mat';
save(FileName,'rho','GS', 'Newton_richardson',...
              'Newton_chebyshev', 'Newton_GS',... 
              'CG', 'richardson','bMMSE', 'chebyshev', 'LAS_bMMSE', 'LAS_Newton_richardson', ...
              'LAS_Newton_chebyshev', 'LAS_Newton_GS', 'LAS_GS', 'LAS_richardson', 'LAS_chebyshev');          
           
          
          
          
          
          
          
fullname = 'results/Plot5_BER_32x128_16QAM_RunID_0.mat';

tmp = load(fullname);
rho                   = tmp.Results.rho';
GS                    = tmp.Results.BER(1,:)';
Newton_richardson     = tmp.Results.BER(2,:)';
Newton_chebyshev      = tmp.Results.BER(3,:)';
Newton_GS             = tmp.Results.BER(4,:)';
CG                    = tmp.Results.BER(5,:)';
richardson            = tmp.Results.BER(6,:)';
bMMSE                 = tmp.Results.BER(7,:)';
chebyshev             = tmp.Results.BER(8,:)';
LAS_bMMSE             = tmp.Results.BER(9,:)';
LAS_Newton_richardson = tmp.Results.BER(10,:)';
LAS_Newton_chebyshev  = tmp.Results.BER(11,:)';
LAS_Newton_GS         = tmp.Results.BER(12,:)';
LAS_GS                = tmp.Results.BER(13,:)';
LAS_richardson        = tmp.Results.BER(14,:)';
LAS_chebyshev         = tmp.Results.BER(15,:)';


FileName = 'Compiled_Data_Origin/plot5_correlation_16QAM_32x128.mat';
save(FileName,'rho','GS', 'Newton_richardson',...
              'Newton_chebyshev', 'Newton_GS',... 
              'CG', 'richardson','bMMSE', 'chebyshev', 'LAS_bMMSE', 'LAS_Newton_richardson', ...
              'LAS_Newton_chebyshev', 'LAS_Newton_GS', 'LAS_GS', 'LAS_richardson', 'LAS_chebyshev');          
          
% \l(1) %(1)
% \l(2) %(2)
% \l(3) %(3)
% \l(4) %(4)
% \l(5) %(5)
% \l(6) %(6)
% \l(7) %(7)
% \l(8) %(8)
% \l(9) %(9)
% \l(10) %(10)
% \l(11) %(11)
% \l(12) %(12)
% \l(13) %(13)
% \l(14) %(14)
% \l(15) %(15)