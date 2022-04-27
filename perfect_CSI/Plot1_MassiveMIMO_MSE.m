% plot 1 Mean square error vs NR/NT ratio plot for different iterations
clc
clear
addpath('results', 'detector', 'sim_profile') ;
rng('default')
par.mod = '64QAM'; % modulation type: 'BPSK','QPSK','16QAM','64QAM'
par.trials = 1000;    par.SNRdB_list = 20;
runId_max = 3;

par.detector = {'Newton_richardson', 'Newton_chebyshev','Newton_GS',...
                 'Conjugate-Gradient', 'richardson',  'chebyshev', 'Gauss-Seidel1'};
par.MT = 32;

par.simName = ['Plot1_MSE_',num2str(par.MT), '_',num2str(par.mod), '_RunID_'];
disp(par.simName)

% parfor i = 0:runId_max
%     Massive_MIMO_all_MSE(par,i)
% end
Massive_MIMO_all_MSE(par,0)




