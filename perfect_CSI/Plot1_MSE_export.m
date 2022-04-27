
%% only MMSE section
clear
fullname = 'results/Plot1_MSE_32_64QAM_RunID_0.mat';
par.detector = {'Newton_richardson', 'Newton_chebyshev','Newton_GS',...
                 'Conjugate-Gradient', 'richardson',  'chebyshev', 'Gauss-Seidel1'};
tmp = load(fullname);
eta_list                 = tmp.Results.eta_temp';
Newton_richardson64QAM   = tmp.Results.MSE(1,:,1)';
Newton_chebyshev64QAM    = tmp.Results.MSE(2,:,1)';
Newton_GS64QAM           = tmp.Results.MSE(3,:,1)';
CG64QAM                  = tmp.Results.MSE(4,:,1)';
richardson64QAM          = tmp.Results.MSE(5,:,1)';
chebyshev64QAM           = tmp.Results.MSE(6,:,1)';
GS64QAM                  = tmp.Results.MSE(7,:,1)';
 
FileName = 'Compiled_Data_Origin/plot1_MSE2_64QAM_32x128.mat';
save(FileName,'eta_list','GS64QAM', 'Newton_richardson64QAM',...
              'Newton_chebyshev64QAM', 'Newton_GS64QAM',... 
              'CG64QAM', 'richardson64QAM','chebyshev64QAM');

          
          
Newton_richardson64QAM   = tmp.Results.MSE(1,:,2)';
Newton_chebyshev64QAM    = tmp.Results.MSE(2,:,2)';
Newton_GS64QAM           = tmp.Results.MSE(3,:,2)';
CG64QAM                  = tmp.Results.MSE(4,:,2)';
richardson64QAM          = tmp.Results.MSE(5,:,2)';
chebyshev64QAM           = tmp.Results.MSE(6,:,2)';
GS64QAM                  = tmp.Results.MSE(7,:,2)';
 
FileName = 'Compiled_Data_Origin/plot1_MSE3_64QAM_32x128.mat';
save(FileName,'eta_list','GS64QAM', 'Newton_richardson64QAM',...
              'Newton_chebyshev64QAM', 'Newton_GS64QAM',... 
              'CG64QAM', 'richardson64QAM','chebyshev64QAM');
          

Newton_richardson64QAM   = tmp.Results.MSE(1,:,3)';
Newton_chebyshev64QAM    = tmp.Results.MSE(2,:,3)';
Newton_GS64QAM           = tmp.Results.MSE(3,:,3)';
CG64QAM                  = tmp.Results.MSE(4,:,3)';
richardson64QAM          = tmp.Results.MSE(5,:,3)';
chebyshev64QAM           = tmp.Results.MSE(6,:,3)';
GS64QAM                  = tmp.Results.MSE(7,:,3)';
 
FileName = 'Compiled_Data_Origin/plot1_MSE4_64QAM_32x128.mat';
save(FileName,'eta_list','GS64QAM', 'Newton_richardson64QAM',...
              'Newton_chebyshev64QAM', 'Newton_GS64QAM',... 
              'CG64QAM', 'richardson64QAM','chebyshev64QAM');
          
