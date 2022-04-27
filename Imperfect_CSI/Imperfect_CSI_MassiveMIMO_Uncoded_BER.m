% uncoded MIMO performance 16 QAM and 64QAM with LAS and other MMSE
% detectors

clc
clear
addpath('results', 'detector', 'sim_profile') ;
% define common simulation parameters
par.trials = 10;    par.block_size = 100;


par.detector = {'Gauss-Seidel1','Newton_richardson', 'Newton_chebyshev','Newton_GS',...
    'Conjugate-Gradient', 'richardson', 'MMSE', 'chebyshev',  ...
    'LAS-bMMSE', 'LAS-Newton_richardson', 'LAS-Newton_chebyshev', 'LAS-Newton_GS',...
    'LAS-GS', 'LAS-richardson','LAS-chebyshev'};
% par.detector = {'Gauss-Seidel1','Newton_richardson', 'Newton_chebyshev','Newton_GS',...
%                'Conjugate-Gradient', 'richardson', 'bMMSE', 'chebyshev','steepest_jacobi','MMSEapprox'};
% %                   'LAS-bMMSE', 'LAS-Newton_richardson', 'LAS-Newton_chebyshev', 'LAS-Newton_GS',...
%                   'LAS-GS', 'LAS-richardson','LAS-chebyshev',};
% par.detector = {'MMSE'};
par.MR = 128; par.MT = 32;     par.PilotNumber = 2*par.MT; par.alg.maxiter = 4;
par.mod = '64QAM'; % modulation type: 'BPSK','QPSK','16QAM','64QAM'
par.SNRdB_list = 12:2:32;
par.band_approx =1;   % 1 means band matrix approximation, 0 --> no approximation
par.w=24;


%% Generate orthogonal pilot symbols
[par.P, par.Pnorm2] = generate_pilot(par.MT,par.PilotNumber);
%% case 1: No correlation
%% Generate correlation matrix
par.correlation = 1; par.zeta = 0; par.theta = pi/3;
[par.Rtx, par.Rrx] = generate_correlation_matrix_2013(par);
par.Rtx = eye(par.MT);
par.simName = ['estimate_BER_No_correlation',num2str(par.MT),'x', num2str(par.MR), '_',num2str(par.mod), '_RunID_'];
disp(par.simName)
Massive_MIMO_all_imperfect_CSI(par,0)

%% case 2: moderate correlation
%% Generate correlation matrix
par.correlation = 1; par.zeta = 0.3; par.theta = pi/3;
[par.Rtx, par.Rrx] = generate_correlation_matrix_2013(par);
par.Rtx = eye(par.MT);
par.simName = ['estimate_BER_moderate_correlation',num2str(par.MT),'x', num2str(par.MR), '_',num2str(par.mod), '_RunID_'];
disp(par.simName)
Massive_MIMO_all_imperfect_CSI(par,0)


























%% 16QAM
%  par.mod = '16QAM'; % modulation type: 'BPSK','QPSK','16QAM','64QAM'
% par.SNRdB_list = 2:2:16;
% par.w=8;
%
% par.simName = ['Plot4_BER_',num2str(par.MT),'x', num2str(par.MR), '_',num2str(par.mod), '_RunID_'];
% disp(par.simName)
%
% % parfor i = 0:runId_max
% %     Massive_MIMO_all(par,i)
% % end
%  Massive_MIMO_all(par,0)
%
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