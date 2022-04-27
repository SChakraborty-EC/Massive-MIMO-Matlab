% effect of omega on LAS performance

clc
clear
addpath('results', 'detector') ;
% define common simulation parameters
par.trials = 1000;    
runId_max = 3;
par.detector = {'LAS-bMMSE', 'LAS-Newton_richardson', 'LAS-GS', ...
                 'LAS-Newton_chebyshev', 'LAS-CG'};

par.MR = 128; par.MT = 32;




% 64QAM modulation
par.SNRdB_list = 19;
par.mod = '64QAM'; % modulation type: 'BPSK','QPSK','16QAM','64QAM'
par.alg.maxiter = 3;       

band_approx =1;   % 1 means band matrix approximation, 0 --> no approximation
par.w=24;

par.simName = ['Plot6_BER_W_',num2str(par.MT),'x', num2str(par.MR), '_',num2str(par.mod), '_RunID_'];
disp(par.simName)

% parfor i = 0:runId_max
%     Massive_MIMO_all_w(par,i)
% end
 Massive_MIMO_all_w(par,0)
 
% 16QAM modulation
 
par.SNRdB_list = 13;
par.mod = '16QAM'; % modulation type: 'BPSK','QPSK','16QAM','64QAM'
par.alg.maxiter = 3;       par.w=24;

par.simName = ['Plot6_BER_W_',num2str(par.MT),'x', num2str(par.MR), '_',num2str(par.mod), '_RunID_'];
disp(par.simName)

% parfor i = 0:runId_max
%     Massive_MIMO_all_w(par,i)
% end
 Massive_MIMO_all_w(par,0)
