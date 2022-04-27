% uncoded MIMO performance 16 QAM and 64QAM with LAS and other MMSE
% detectors

clc
clear
addpath('results', 'detector') ;
% define common simulation parameters
par.trials = 80000;   
runId_max = 3;

% par.detector = {'Gauss-Seidel1','Newton_richardson', 'Newton_chebyshev','Newton_GS',...
%                  'Conjugate-Gradient', 'richardson', 'bMMSE', 'chebyshev',  ...
%                   'LAS-bMMSE', 'LAS-Newton_richardson', 'LAS-Newton_chebyshev', 'LAS-Newton_GS',... 
%                   'LAS-GS', 'LAS-richardson','LAS-chebyshev',};
par.detector = {'Gauss-Seidel1','Newton_richardson', 'Newton_chebyshev','Newton_GS',...
               'Conjugate-Gradient', 'richardson', 'bMMSE', 'chebyshev','steepest_jacobi','MMSEapprox'};
%                   'LAS-bMMSE', 'LAS-Newton_richardson', 'LAS-Newton_chebyshev', 'LAS-Newton_GS',... 
%                   'LAS-GS', 'LAS-richardson','LAS-chebyshev',};

par.MR = 128; par.MT = 32;
par.alg.maxiter = 3;

%% 64QAM
par.mod = '64QAM'; % modulation type: 'BPSK','QPSK','16QAM','64QAM'
par.SNRdB_list = 6:2:28;
par.band_approx =1;   % 1 means band matrix approximation, 0 --> no approximation
par.w=16;

par.simName = ['Plot4_BER_',num2str(par.MT),'x', num2str(par.MR), '_',num2str(par.mod), '_RunID_'];
disp(par.simName)
Massive_MIMO_all(par,0)

%% 16QAM 
 par.mod = '16QAM'; % modulation type: 'BPSK','QPSK','16QAM','64QAM'
par.SNRdB_list = 2:2:18;
par.w=8;

par.simName = ['Plot4_BER_',num2str(par.MT),'x', num2str(par.MR), '_',num2str(par.mod), '_RunID_'];
disp(par.simName)

Massive_MIMO_all(par,0)
% --------------------------------------------------------------
%% Antenna 64 --------------------------------------------------
% --------------------------------------------------------------
par.MR = 256; par.MT = 64;
par.alg.maxiter = 3;

%% 64QAM
par.mod = '64QAM'; % modulation type: 'BPSK','QPSK','16QAM','64QAM'
par.SNRdB_list = 6:2:28;
par.band_approx =1;   % 1 means band matrix approximation, 0 --> no approximation
par.w=16;

par.simName = ['Plot4_BER_',num2str(par.MT),'x', num2str(par.MR), '_',num2str(par.mod), '_RunID_'];
disp(par.simName)
Massive_MIMO_all(par,0)

%% 16QAM 
 par.mod = '16QAM'; % modulation type: 'BPSK','QPSK','16QAM','64QAM'
par.SNRdB_list = 2:2:18;
par.w=8;

par.simName = ['Plot4_BER_',num2str(par.MT),'x', num2str(par.MR), '_',num2str(par.mod), '_RunID_'];
disp(par.simName)

Massive_MIMO_all(par,0)
%  clear
%  fullname = ['results/Plot4_BER_16x64_64QAM_RunID_0.mat'];
%  
%  tmp = load(fullname);
%  tmp.Results.BER
%  close all
% marker_style = {'bo-','rs--','mv-.','kp:','g*-','c>--','yx:', 'kx-',...
%                 'kh:', 'b+--', 'rs-', 'gd:', 'ys--', 'c+--', 'k*', 'kh:'};
% figure(1)
% for d=1:length(tmp.Results.detector)
%     if d==1
%         semilogy(tmp.Results.SNRdB_list,tmp.Results.BER(d,:),marker_style{d},'LineWidth',2)
%         hold on
%     else
%         
%         semilogy(tmp.Results.SNRdB_list,tmp.Results.BER(d,:),marker_style{d},'LineWidth',2)
%     end
% end
% hold off
% grid on
% xlabel('average SNR per receive antenna [dB]','FontSize',12)
% ylabel('bit error rate (BER)','FontSize',12)
% axis([min(tmp.Results.SNRdB_list) max(tmp.Results.SNRdB_list) 1e-6 1])
% legend(tmp.Results.detector,'FontSize',12)
% set(gca,'FontSize',12)
% 
%  fullname = ['results/Plot4_BER_16x64_16QAM_RunID_0.mat'];
%  
%  tmp = load(fullname);
%  tmp.Results.BER
% 
% marker_style = {'bo-','rs--','mv-.','kp:','g*-','c>--','yx:', 'kx-',...
%                 'kh:', 'b+--', 'rs-', 'gd:', 'ys--', 'c+--', 'k*', 'kh:'};
% figure(2)
% for d=1:length(tmp.Results.detector)
%     if d==1
%         semilogy(tmp.Results.SNRdB_list,tmp.Results.BER(d,:),marker_style{d},'LineWidth',2)
%         hold on
%     else
%         
%         semilogy(tmp.Results.SNRdB_list,tmp.Results.BER(d,:),marker_style{d},'LineWidth',2)
%     end
% end
% hold off
% grid on
% xlabel('average SNR per receive antenna [dB]','FontSize',12)
% ylabel('bit error rate (BER)','FontSize',12)
% axis([min(tmp.Results.SNRdB_list) max(tmp.Results.SNRdB_list) 1e-6 1])
% legend(tmp.Results.detector,'FontSize',12)
% set(gca,'FontSize',12)