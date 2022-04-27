% Effect of spetial correletion 
% for fixed SNR BER vs Spetial correlation coefficient plot

clc
clear
addpath('results', 'detector') ;
% define common simulation parameters
par.trials = 1000;    
runId_max = 2;

par.detector = {'Gauss-Seidel1','Newton_richardson', 'Newton_chebyshev','Newton_GS',...
                 'Conjugate-Gradient', 'richardson', 'bMMSE', 'chebyshev', ...
                'LAS-bMMSE', 'LAS-Newton_richardson', 'LAS-Newton_chebyshev', 'LAS-Newton_GS',... 
                'LAS-GS', 'LAS-richardson','LAS-chebyshev'};             
par.MR = 128; par.MT = 32;
par.alg.maxiter = 3;
par.rho_max = 0.5;

% 64QAM modulation
par.mod = '64QAM'; % modulation type: 'BPSK','QPSK','16QAM','64QAM'
par.band_approx =1;   % 1 means band matrix approximation, 0 --> no approximation
par.w=20;

par.SNRdB_list = 18;
par.simName = ['Plot5_BER_W20_',num2str(par.MT),'x', num2str(par.MR), '_',num2str(par.mod), '_RunID_'];
disp(par.simName)

% parfor i = 0:runId_max
%     Massive_MIMO_all_correlated_variation(par,i)
% end
 Massive_MIMO_all_correlated_variation(par,0)

 
 
% 16QAM modulation 
par.mod = '16QAM'; % modulation type: 'BPSK','QPSK','16QAM','64QAM'
par.band_approx =0;   % 1 means band matrix approximation, 0 --> no approximation
par.w=4;
par.SNRdB_list = 13;
par.simName = ['Plot5_BER_',num2str(par.MT),'x', num2str(par.MR), '_',num2str(par.mod), '_RunID_'];
disp(par.simName)

% parfor i = 0:runId_max
%     Massive_MIMO_all_correlated_variation(par,i)
% end
%  Massive_MIMO_all_correlated_variation(par,0)
