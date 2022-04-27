% BER performance for the variation of base station antennas for different
% algorithms

clc
clear
addpath('results', 'detector') ;
% define common simulation parameters
par.trials = 1000;    
runId_max = 2;

par.detector = {'Gauss-Seidel1','Newton_richardson', 'Newton_chebyshev','Newton_GS',...
                 'Conjugate-Gradient', 'richardson', 'bMMSE', 'chebyshev'};             
par.MR = 128; 
par.alg.maxiter = 4;


% 64QAM modulation
par.mod = '64QAM'; % modulation type: 'BPSK','QPSK','16QAM','64QAM'

par.SNRdB_list = 16;
par.simName = ['Plot7_BER_tx_variation16dB4itr_',num2str(par.MR), '_',num2str(par.mod), '_RunID_'];
disp(par.simName)

% parfor i = 0:runId_max
%     Massive_MIMO_all_correlated_variation(par,i)
% end
 Massive_MIMO_all_tx_variation(par,0)

 
