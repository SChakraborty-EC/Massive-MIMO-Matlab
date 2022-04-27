% -----------------------------------------------------
% -----------------------------------------------------

function Massive_MIMO_all_MSE(par, runId)
% use runId random seed (enables reproducibility)
par.runId=runId;
rng(par.runId);
eta_max = 10;
% set up Gray-mapped constellation alphabet (according to IEEE 802.11)
switch (par.mod)
    case 'BPSK',
        par.symbols = [ -1 1 ];
        
    case 'QPSK',
        par.symbols = [ -1-1i,-1+1i, ...
            +1-1i,+1+1i ];
        par.RealSymbols = [-1 1];
        
    case '16QAM',
        par.symbols = [ -3-3i,-3-1i,-3+3i,-3+1i, ...
            -1-3i,-1-1i,-1+3i,-1+1i, ...
            +3-3i,+3-1i,+3+3i,+3+1i, ...
            +1-3i,+1-1i,+1+3i,+1+1i ];
        
        par.M = 3;
        
        par.RealSymbols = [-3 -1 3 1];
    case '64QAM',
        par.symbols = [ -7-7i,-7-5i,-7-1i,-7-3i,-7+7i,-7+5i,-7+1i,-7+3i, ...
            -5-7i,-5-5i,-5-1i,-5-3i,-5+7i,-5+5i,-5+1i,-5+3i, ...
            -1-7i,-1-5i,-1-1i,-1-3i,-1+7i,-1+5i,-1+1i,-1+3i, ...
            -3-7i,-3-5i,-3-1i,-3-3i,-3+7i,-3+5i,-3+1i,-3+3i, ...
            +7-7i,+7-5i,+7-1i,+7-3i,+7+7i,+7+5i,+7+1i,+7+3i, ...
            +5-7i,+5-5i,+5-1i,+5-3i,+5+7i,+5+5i,+5+1i,+5+3i, ...
            +1-7i,+1-5i,+1-1i,+1-3i,+1+7i,+1+5i,+1+1i,+1+3i, ...
            +3-7i,+3-5i,+3-1i,+3-3i,+3+7i,+3+5i,+3+1i,+3+3i ];
        
        par.RealSymbols = [-7 -5 -1 -3 7 5 1 3];
        par.M = 7;
        
    case '256QAM',
        par.RealSymbols = [-15 -13 -9 -11 -1 -3 -7 -5 15 13 9 11 1 3 7 5];
        par.symbols = zeros(1,256);
        index = 1;
        for i = 1:16
            for j = 1:16
                par.symbols (1,index)  = par.RealSymbols (i) +1i*par.RealSymbols (j);
                index = index+1;
            end
        end
        par.M =  15;
        
end

% extract average symbol energy
par.Es = mean(abs(par.symbols).^2);

% precompute bit labels
par.Q = log2(length(par.symbols)); % number of bits per symbol
par.bits = de2bi(0:length(par.symbols)-1,par.Q,'left-msb');

% precompute real bit labels
par.realQ = log2(length(par.RealSymbols)); % number of bits per symbol
par.realbits = de2bi(0:length(par.RealSymbols)-1,par.realQ,'left-msb');

% track simulation time
time_elapsed = 0;

% -- start simulation

% initialize result arrays (detector x SNR)
res.MSE = zeros(length(par.detector),eta_max,3); % bit error rate

% generate random bit stream (antenna x bit x trial)
bits = randi([0 1],par.MT,par.Q,par.trials);

% trials loop
tic
for eta = 1:eta_max
    par.MR = eta*par.MT;
for t=1:par.trials
    
    % generate transmit symbol
    idx = bi2de(bits(:,:,t),'left-msb')+1;
    s = par.symbols(idx).';
    
    % generate iid Gaussian channel matrix & noise vector
    n = sqrt(0.5)*(randn(par.MR,1)+1i*randn(par.MR,1));
    H = sqrt(0.5)*(randn(par.MR,par.MT)+1i*randn(par.MR,par.MT));
    
    % transmit over noiseless channel (will be used later)
    x = H*s;
    % compute noise variance (average SNR per receive antenna is: SNR=MT*Es/N0)
    N0 = par.MT*par.Es*10^(-par.SNRdB_list/10);
    
    % transmit data over noisy channel
    y = x+sqrt(N0)*n;
    xhat_MMSE = bMMSE(par,H,y,N0);
    % SNR loop
    for k=1:3
        par.alg.maxiter = k+1;       
        for d=1:length(par.detector)
            % algorithm loop
            switch (par.detector{d}) % select algorithms
                case 'Gauss-Seidel1' % Gauss-Seidel detector
                    xhat = Gauss_Seidel1(par,H,y,N0);

                case 'Neumann' % coordinate descent
                    xhat = Neumann(par,H,y,N0);

                case 'Conjugate-Gradient' % conjugate gradient detector
                    xhat = CG(par,H,y,N0);

                case 'richardson' % conjugate gradient detector
                    xhat = richardson(par,H,y,N0);

                case 'SOR' % conjugate gradient detector
                    xhat = SOR(par,H,y,N0);

                case 'Newton_iteration' % conjugate gradient detector
                    xhat = Newton_iteration(par,H,y,N0);

                case 'jacobi' % conjugate gradient detector
                    xhat = jacobi(par,H,y,N0);

                case 'chebyshev' % conjugate gradient detector
                    xhat = chebyshev(par,H,y,N0);

                case 'Newton_richardson' % conjugate gradient detector
                    xhat = Newton_richardson(par,H,y,N0);
                    
                case 'GS_richardson' % conjugate gradient detector
                    xhat = GS_richardson(par,H,y,N0);
                    
                case 'Newton_chebyshev' % conjugate gradient detector
                    xhat = Newton_chebyshev(par,H,y,N0);
                case 'Newton_SOR' % conjugate gradient detector
                    xhat = Newton_SOR(par,H,y,N0);
                case 'MMSEapprox', % biased MMSE detector
                    xhat = MMSE_approx(par,H,y, N0);
                case 'Newton_GS' % conjugate gradient detector
                    xhat = Newton_GS(par,H,y,N0);

                otherwise,
                    error('par.detector type not defined.')
            end

            res.MSE(d,eta,k) = res.MSE(d,eta,k) + norm(xhat - xhat_MMSE)^2;
        end % detector loop
    end % SNR loop
    
    % keep track of simulation time
    if toc>10
        time=toc;
        time_elapsed = time_elapsed + time;
        fprintf('estimated remaining simulation time: %3.0f min.\n',time_elapsed*(par.trials/t-1)/60);
        tic
    end
    %                                 fprintf('packet =%0.5g,Symbol =%0.5g, SNR (dB) = %0.5g , Error = %0.5g \n',pack,symb,SNR_dB,tmp);
    
end % trials loop
end
% -- simulation completed!
% disp(' =====================================================');
% disp(sprintf(' ==================== Total Simulation Time: %f sec. ',Result.total_time));
% disp(' =====================================================');
% -- show results (generates fairly nice Matlab plot)

FileName = sprintf('Results/%s%d.mat',par.simName, par.runId);
if exist(FileName, 'file')
    delete(FileName);
    fprintf('\n %s%d.mat  file is deleted\n', par.simName,  par.runId)
end

% -- save final results (par and res structure)

% save([ par.simName '_' num2str(par.runId) ],'par','res');
Results.MSE=res.MSE/(par.trials*par.MT^2);
Results.time_elapsed = time_elapsed;
Results.simname= par.simName;
Results.SNRdB_list = par.SNRdB_list ;
Results.trials =  par.trials;
Results.Detector_list  = par.detector;
Results.FileName = sprintf('results/%s%d.mat',par.simName, par.runId);

Results.MT = par.MT;
Results.MR = par.MT*[1:eta_max];
Results.eta_temp = [1:eta_max];
save(Results.FileName,'Results');
eta_temp = [1:eta_max];
close all
marker_style = {'bo-','rs--','mv-.','kp:','g*-','c>--','yx:', 'kx-',...
    'kh:', 'b+--', 'rs-', 'gd:', 'ys--', 'c+--', 'k*', 'kh:'};
figure(1)
for d=1:length(par.detector)
    if d==1
        semilogy(eta_temp,Results.MSE(d,:,1),marker_style{d},'LineWidth',2)
        hold on
    else
        
        semilogy(eta_temp,Results.MSE(d,:,1),marker_style{d},'LineWidth',2)
    end
end

hold off
grid on
xlabel('number of receive antenna [dB]','FontSize',12)
ylabel('MSE','FontSize',12)
% axis([min(Results.MR) max(Results.MR) 0 1e-4])
legend(par.detector,'FontSize',12)
set(gca,'FontSize',12)


figure(2)
for d=1:length(par.detector)
    if d==1
        semilogy(eta_temp,Results.MSE(d,:,2),marker_style{d},'LineWidth',2)
        hold on
    else
        
        semilogy(eta_temp,Results.MSE(d,:,2),marker_style{d},'LineWidth',2)
    end
end

hold off
grid on
xlabel('number of receive antenna [dB]','FontSize',12)
ylabel('MSE','FontSize',12)
% axis([min(Results.MR) max(Results.MR) 0 1e-4])
legend(par.detector,'FontSize',12)
set(gca,'FontSize',12)

figure(3)
for d=1:length(par.detector)
    if d==1
        semilogy(eta_temp,Results.MSE(d,:,3),marker_style{d},'LineWidth',2)
        hold on
    else
        
        semilogy(eta_temp,Results.MSE(d,:,3),marker_style{d},'LineWidth',2)
    end
end

hold off
grid on
xlabel('number of receive antenna [dB]','FontSize',12)
ylabel('MSE','FontSize',12)
% axis([min(Results.MR) max(Results.MR) 0 1e-4])
legend(par.detector,'FontSize',12)
set(gca,'FontSize',12)
end
