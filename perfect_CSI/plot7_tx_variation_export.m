clear
fullname = 'results/Plot7_BER_tx_variation16dB_128_64QAM_RunID_0.mat';
par.detector = {'Gauss-Seidel1','Newton_richardson', 'Newton_chebyshev','Newton_GS',...
                 'Conjugate-Gradient', 'richardson', 'bMMSE', 'chebyshev'};             
tmp = load(fullname);
tx                      = tmp.Results.tx';
GS64QAM                  = tmp.Results.BER(1,:)';
Newton_richardson64QAM   = tmp.Results.BER(2,:)';
Newton_chebyshev64QAM    = tmp.Results.BER(3,:)';
Newton_GS64QAM           = tmp.Results.BER(4,:)';
CG64QAM                  = tmp.Results.BER(5,:)';
richardson64QAM          = tmp.Results.BER(6,:)';
bMMSE                    = tmp.Results.BER(7,:)';
chebyshev64QAM           = tmp.Results.BER(8,:)';
 
 
 
 
FileName = 'Compiled_Data_Origin/plot7_64QAM_32x128.mat';
save(FileName,'tx','GS64QAM', 'Newton_richardson64QAM',...
              'Newton_chebyshev64QAM', 'Newton_GS64QAM',... 
              'CG64QAM', 'richardson64QAM','bMMSE','chebyshev64QAM');
          
          

              
