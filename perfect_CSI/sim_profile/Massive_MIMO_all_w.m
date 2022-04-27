% -----------------------------------------------------
% -----------------------------------------------------

function Massive_MIMO_all_w(par, runId)
% use runId random seed (enables reproducibility)
par.runId=runId;
rng(par.runId);
w_max = 2*par.MT-8;
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
res.BER = zeros(length(par.detector),w_max,3); % bit error rate

% generate random bit stream (antenna x bit x trial)
bits = randi([0 1],par.MT,par.Q,par.trials);

% trials loop
tic
for w = 1:2:w_max
    par.w=w-1;
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
    % SNR loop
    for k=1:3
        par.alg.maxiter = k+1;       
        for d=1:length(par.detector)
            % algorithm loop
            switch (par.detector{d}) % select algorithms
                case 'LAS-bMMSE', % ML detection using sphere decoding
                    par.MMSE_init = 'bMMSE';
                    bithat = LAS1(par,H,y, N0);
                    
                case 'LAS-CG', % ML detection using sphere decoding
                    par.MMSE_init = 'Conjugate-Gradient';
                    bithat = LAS1(par,H,y, N0);
                case 'LAS-GS', % ML detection using sphere decoding
                    par.MMSE_init = 'Gauss-Seidel';
                    bithat = LAS1(par,H,y, N0);

                case 'LAS-richardson', % ML detection using sphere decoding
                    par.MMSE_init = 'richardson';
                    bithat = LAS1(par,H,y, N0);
                case 'LAS-SOR', % ML detection using sphere decoding
                    par.MMSE_init = 'SOR';
                    bithat = LAS1(par,H,y, N0);
                case 'LAS-Newton_iteration', % ML detection using sphere decoding
                    par.MMSE_init = 'Newton_iteration';
                    bithat = LAS1(par,H,y, N0);
                case 'LAS-jacobi', % ML detection using sphere decoding
                    par.MMSE_init = 'jacobi';
                    bithat = LAS1(par,H,y, N0);
                case 'LAS-Newton_richardson', % ML detection using sphere decoding
                    par.MMSE_init = 'Newton_richardson';
                    bithat = LAS1(par,H,y, N0);
                case 'LAS-Newton_chebyshev', % ML detection using sphere decoding
                    par.MMSE_init = 'Newton_chebyshev';
                    bithat = LAS1(par,H,y, N0);
                otherwise,
                    error('par.detector type not defined.')
            end
            res.BER(d,w,k) = res.BER(d,w,k) + sum(sum(bits(:,:,t)~=bithat))/(par.MT*par.Q);
%             res.MSE(d,eta,k) = res.MSE(d,eta,k) + norm(xhat - xhat_MMSE)^2;
        end % detector loop
    end % SNR loop
    
%     % keep track of simulation time
%     if toc>10
%         time=toc;
%         time_elapsed = time_elapsed + time;
%         fprintf('estimated remaining simulation time: %3.0f min.\n',time_elapsed*(par.trials/t-1)/60);
%         tic
%     end
    %                                 fprintf('packet =%0.5g,Symbol =%0.5g, SNR (dB) = %0.5g , Error = %0.5g \n',pack,symb,SNR_dB,tmp);
    
end % trials loop

    % keep track of simulation time
    if toc>10
        time=toc;
        time_elapsed = time_elapsed + time;
        fprintf('estimated remaining simulation time: %3.0f min.\n',time_elapsed*(w_max/w-1)/60);
        tic
    end

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
Results.BER=res.BER/par.trials;
Results.time_elapsed = time_elapsed;
Results.simname= par.simName;
Results.SNRdB_list = par.SNRdB_list ;
Results.trials =  par.trials;

Results.FileName = sprintf('results/%s%d.mat',par.simName, par.runId);
save(Results.FileName,'Results');

Results.MR = par.MR;
w_temp = [0:2:w_max-1];
close all
marker_style = {'bo-','rs--','mv-.','kp:','g*-','c>--','yx:', 'kx-',...
    'kh:', 'b+--', 'rs-', 'gd:', 'ys--', 'c+--', 'k*', 'kh:'};
figure(1)
for d=1:length(par.detector)
    if d==1
        semilogy(w_temp,Results.BER(d,:,1),marker_style{d},'LineWidth',2)
        hold on
    else
        
        semilogy(w_temp,Results.BER(d,:,1),marker_style{d},'LineWidth',2)
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
        semilogy(w_temp,Results.BER(d,:,2),marker_style{d},'LineWidth',2)
        hold on
    else
        
        semilogy(w_temp,Results.BER(d,:,2),marker_style{d},'LineWidth',2)
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
        semilogy(w_temp,Results.BER(d,:,3),marker_style{d},'LineWidth',2)
        hold on
    else
        
        semilogy(w_temp,Results.BER(d,:,3),marker_style{d},'LineWidth',2)
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
