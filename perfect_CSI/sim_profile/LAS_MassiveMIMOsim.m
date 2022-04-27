% -----------------------------------------------------
% -----------------------------------------------------

function LAS_MassiveMIMOsim(par,runId)

% -- set up default/custom parameters

   
    % set default simulation parameters
%     par.simName = 'ERR_4x4_16QAM'; % simulation name (used for saving results)
     par.runId = runId; % simulation ID (used to reproduce results)
%     par.MR = 256; % receive antennas
%     par.MT = 64; % transmit antennas (set not larger than MR!)
%     par.mod = '16QAM'; % modulation type: 'BPSK','QPSK','16QAM','64QAM'
%     par.trials = 1000; % number of Monte-Carlo trials (transmissions)
%     par.SNRdB_list = 5:2:15; % list of SNR [dB] values to be simulated
%     par.detector = {'Gauss-Seidel','bMMSE', 'MMSEapprox', 'MMSEapprox2','MMSEapprox2_real'}; % define detector(s) to be simulated

%         par.detector = {'Gauss-Seidel','bMMSE', 'MMSEapprox', 'MMSEapprox2','MMSE_approx_exact_eigen', 'Neumann','OCDBOX','ADMIN', 'Conjugate-Gradient'}; % define detector(s) to be simulated
%     par.detector = {'LAS1', 'bMMSE', 'Gauss-Seidel'}; % define detector(s) to be simulated
%      par.detector ={'TABU-CG','TABU-NU','TABU-GS','TABU-bMMSE','LAS-CG','LAS-NU','LAS-GS','LAS-bMMSE','Conjugate-Gradient','Neumann','Gauss-Seidel','OCDBOX','ADMIN', 'MMSEapprox', 'bMMSE'};
%     par.detector ={'LAS-proposed','LAS-CG','LAS-NU','LAS-GS','LAS-bMMSE','Conjugate-Gradient','Neumann','Gauss-Seidel','OCDBOX','ADMIN', 'MMSEapprox2_real', 'bMMSE'};
    
    par.alg.maxiter = 3;

% -- initialization

% use runId random seed (enables reproducibility)
rng(par.runId);

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
        
        par.RealSymbols = [-3 -1 3 1];
        par.M =  3;
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
        par.M =  7;
        
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
res.BER = zeros(length(par.detector),length(par.SNRdB_list)); % bit error rate

% generate random bit stream (antenna x bit x trial)
bits = randi([0 1],par.MT,par.Q,par.trials);

% trials loop
tic
for t=1:par.trials
    
    % generate transmit symbol
    idx = bi2de(bits(:,:,t),'left-msb')+1;
    s = par.symbols(idx).';
    
    % generate iid Gaussian channel matrix & noise vector
    n = sqrt(0.5)*(randn(par.MR,1)+1i*randn(par.MR,1));
    H = sqrt(0.5)*(randn(par.MR,par.MT)+1i*randn(par.MR,par.MT));
    
    % transmit over noiseless channel (will be used later)
    x = H*s;
    
    % SNR loop
    for k=1:length(par.SNRdB_list)
        
        % compute noise variance (average SNR per receive antenna is: SNR=MT*Es/N0)
        N0 = par.MT*par.Es*10^(-par.SNRdB_list(k)/10);
        par.N0 = N0;
        % transmit data over noisy channel
        y = x+sqrt(N0)*n;
        %         s1 = [real(s)' imag(s)']'
        
        % algorithm loop
        %         original_norm = norm(y - H*s)^2
        for d=1:length(par.detector)
            
            switch (par.detector{d}) % select algorithms
                case 'ZF', % zero-forcing detection
                    [idxhat,bithat] = ZF(par,H,y);
                case 'bMMSE', % biased MMSE detector
                    [idxhat,bithat] = bMMSE(par,H,y,N0);
                case 'uMMSE', % unbiased MMSE detector
                    [idxhat,bithat] = uMMSE(par,H,y,N0);
                case 'ML', % ML detection using sphere decoding
                    [idxhat,bithat] = ML(par,H,y);
                case 'TABU', % ML detection using sphere decoding
                    [idxhat, bithat] = tabu_search_modified(par,H,y);
                case 'VBLAST', % ML detection using sphere decoding
                    [idxhat, bithat] = VBLAST(par,H,y);
                case 'LAS-CG', % ML detection using sphere decoding
                    par.MMSE_init = 'Conjugate-Gradient';
                    bithat = LAS1(par,H,y);
                case 'LAS-GS', % ML detection using sphere decoding
                    par.MMSE_init = 'Gauss-Seidel';
                    bithat = LAS1(par,H,y);
                case 'LAS-NU', % ML detection using sphere decoding
                    par.MMSE_init = 'Neumann';
                    bithat = LAS1(par,H,y);
                    
                case 'LAS-proposed', % ML detection using sphere decoding
                    par.MMSE_init = 'MMSE_approx2_real';
                    bithat = LAS1(par,H,y);
                case 'LAS-bMMSE', % ML detection using sphere decoding
                    par.MMSE_init = 'bMMSE';
                    bithat = LAS1(par,H,y);
                case 'TABU-CG', % ML detection using sphere decoding
                    par.MMSE_init = 'Conjugate-Gradient';
                    bithat = tabu_search_simplefied_simplified(par,H,y);
                case 'TABU-GS', % ML detection using sphere decoding
                    par.MMSE_init = 'Gauss-Seidel';
                    bithat = tabu_search_simplefied_simplified(par,H,y);
                case 'TABU-NU', % ML detection using sphere decoding
                    par.MMSE_init = 'Neumann';
                    bithat = tabu_search_simplefied_simplified(par,H,y);
                case 'TABU-bMMSE', % ML detection using sphere decoding
                    par.MMSE_init = 'bMMSE';
                    bithat = tabu_search_simplefied_simplified(par,H,y);
                case 'LLAS', % ML detection using sphere decoding
                    [idxhat, bithat] = layered_LAS(par,H,y);
                case 'Gauss-Seidel' % Gauss-Seidel detector
                    [idxhat,bithat] = Gauss_Seidel(par,H,y,N0);
                case 'MF' % Matched Filter
                    [idxhat,bithat] = MF(par,H,y,N0);
                case 'MMSE' % MMSE detector
                    [idxhat,bithat] = MMSE(par,H,y,N0);
                case 'SIMO' % SIMO lower bound
                    [idxhat,bithat] = SIMO(par,H,y,N0,s);
                case 'ADMIN' % ADMM-based Infinity Norm detector
                    [idxhat,bithat] = ADMIN(par,H,y,N0);
                case 'OCDBOX' % co-ordinate descent (optimized) detector
                    [idxhat,bithat] = OCDBOX(par,H,y);
                case 'Neumann' % coordinate descent
                    [idxhat,bithat] = Neumann(par,H,y,N0);
                case 'Gauss-Seidel' % Gauss-Seidel detector
                    [idxhat,bithat] = Gauss_Seidel(par,H,y,N0);
                case 'Conjugate-Gradient' % conjugate gradient detector
                    [idxhat,bithat] = CG(par,H,y,N0);
                case 'MMSEapprox', % biased MMSE detector
                    [idxhat,bithat] = MMSE_approx(par,H,y, N0);
                 case 'MMSEapprox2', % biased MMSE detector
                    [idxhat,bithat] = MMSE_approx2(par,H,y, N0);
                case 'MMSEapprox2_real', % biased MMSE detector
                    [idxhat,bithat] = MMSE_approx2_real(par,H,y, N0);

                 case 'MMSE_approx_exact_eigen', % biased MMSE detector
                    [idxhat,bithat] = MMSE_approx_exact_eigen(par,H,y, N0);
                   
                case 'bMMSE', % biased MMSE detector
                    [idxhat,bithat] = bMMSE(par,H,y,N0);
                    
                otherwise,
                    error('par.detector type not defined.')
            end
            res.BER(d,k) = res.BER(d,k) + sum(sum(bits(:,:,t)~=bithat))/(par.MT*par.Q);
        end % algorithm loop
        
    end % SNR loop
    
    % keep track of simulation time
    if toc>10
        time=toc;
        time_elapsed = time_elapsed + time;
        fprintf('estimated remaining simulation time: %3.0f min.\n',time_elapsed*(par.trials/t-1)/60);
        tic
    end
    
end % trials loop

% normalize results
res.BER = res.BER/par.trials;
res.time_elapsed = time_elapsed;

% -- save final results (par and res structure)

save([ par.simName '_' num2str(par.runId) ],'par','res');
close all
% -- show results (generates fairly nice Matlab plot)
marker_style = {'bo-','rs--','mv-.','kp:','g*-','c>--','yx:', 'kx-', 'kh:', 'b+--', 'rs-', 'gd:', 'ys--', 'c+--', 'k*', 'kh:'};
figure(1)
for d=1:length(par.detector)
    if d==1
        semilogy(par.SNRdB_list,res.BER(d,:),marker_style{d},'LineWidth',2)
        hold on
    else
        semilogy(par.SNRdB_list,res.BER(d,:),marker_style{d},'LineWidth',2)
    end
end

hold off
grid on
xlabel('average SNR per receive antenna [dB]','FontSize',12)
ylabel('bit error rate (BER)','FontSize',12)
axis([min(par.SNRdB_list) max(par.SNRdB_list) 1e-6 1])
legend(par.detector,'FontSize',12)
set(gca,'FontSize',12)
% figure(2)


end

% -- set of detector functions

%% zero-forcing (ZF) detector
function [idxhat,bithat] = ZF(par,H,y)
xhat = H\y;
[~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
bithat = par.bits(idxhat,:);
end

%% biased MMSE detector (bMMSE)
function [idxhat,bithat] = bMMSE(par,H,y,N0)
xhat = (H'*H+(N0/par.Es)*eye(par.MT))\(H'*y);
[~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
bithat = par.bits(idxhat,:);
end

%% unbiased MMSE detector (uMMSE)
function [idxhat,bithat] = uMMSE(par,H,y,N0)
W = (H'*H+(N0/par.Es)*eye(par.MT))\(H');
xhat = W*y;
G = real(diag(W*H));
[~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-G*par.symbols).^2,[],2);
bithat = par.bits(idxhat,:);
end

%% ML detection using sphere decoding
function [idxML,bitML] = ML(par,H,y)

% -- initialization
Radius = inf;
PA = zeros(par.MT,1); % path
ST = zeros(par.MT,length(par.symbols)); % stack

% -- preprocessing
[Q,R] = qr(H,0);
y_hat = Q'*y;

% -- add root node to stack
Level = par.MT;
ST(Level,:) = abs(y_hat(Level)-R(Level,Level)*par.symbols.').^2;

% -- begin sphere decoder
while ( Level<=par.MT )
    % -- find smallest PED in boundary
    [minPED,idx] = min( ST(Level,:) );
    
    % -- only proceed if list is not empty
    if minPED<inf
        ST(Level,idx) = inf; % mark child as tested
        NewPath = [ idx ; PA(Level+1:end,1) ]; % new best path
        
        % -- search child
        if ( minPED<Radius )
            % -- valid candidate found
            if ( Level>1 )
                % -- expand this best node
                PA(Level:end,1) = NewPath;
                Level = Level-1; % downstep
                DF = R(Level,Level+1:end) * par.symbols(PA(Level+1:end,1)).';
                ST(Level,:) = minPED + abs(y_hat(Level)-R(Level,Level)*par.symbols.'-DF).^2;
            else
                % -- valid leaf found
                idxML = NewPath;
                bitML = par.bits(idxML',:);
                % -- update radius (radius reduction)
                Radius = minPED;
            end
        end
    else
        % -- no more childs to be checked
        Level=Level+1;
    end
end

end
% -- set of detector functions

%% Matched filter

function [idxhat,bithat] = MF(par,H,y)

xhat = H' * y / norm(H(:));
[~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
bithat = par.bits(idxhat,:);
end

%% MMSE detector (MMSE)
function [idxhat,bithat] = MMSE(par,H,y,N0)
xhat = (H'*H+(N0/par.Es)*eye(par.MT))\(H'*y);
[~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
bithat = par.bits(idxhat,:);
end

%% SIMO bound
function [idxhat,bithat] = SIMO(par,H,y,s)
z = y-H*s;
xhat = zeros(par.MT,1);
for m=1:par.MT
    hm = H(:,m);
    yhat = z+hm*s(m,1);
    xhat(m,1) = hm'*yhat/norm(hm,2)^2;
end
[~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
bithat = par.bits(idxhat,:);
end


%% Neumann-Series Approximation based massive MIMO detection

function [idxhat,bithat] = Neumann(par,H,y,N0)
A = H'*H+(N0/par.Es)*eye(par.MT);
MF = H'*y;

D = diag(diag(A));
E = triu(A,1)+tril(A,-1);
Ainv = 0;
for i = 0:par.alg.maxiter
    Ainv = Ainv+((-inv(D)*E)^i)*inv(D);
end

xhat = Ainv*MF;
[~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
bithat = par.bits(idxhat,:);
end

%% Gauss-Seidel massive MIMO detection

function [idxhat,bithat] = Gauss_Seidel(par,H,y,N0)
A = H'*H+(N0/par.Es)*eye(par.MT);
MF = H'*y;

D = diag(diag(A));
E = -triu(A,1);
F = -tril(A,-1);

xhat = diag(inv(D));% inv(D)*MF;  %%% Check Gauss Seidel detection paper
for i = 0:par.alg.maxiter
    xhat = inv(D-E)*(F*xhat+MF);
end

[~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
bithat = par.bits(idxhat,:);
end



%% Conjugate Gradient massive MIMO detection

function [idxhat,bithat] = CG(par,H,y,N0)
A = H'*H+(N0/par.Es)*eye(par.MT);
MF = H'*y;

r = MF;
p = r;
v = zeros(par.MT,1);

for k = 1:par.alg.maxiter
    e = A*p;
    alpha  = norm(r)^2/(p'*e);
    v = v+alpha*p;
    new_r = r-alpha*e;
    beta = norm(new_r)^2/norm(r)^2;
    p = new_r+beta*p;
    r = new_r;
end

xhat = v;

[~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
bithat = par.bits(idxhat,:);

end


%% ADMM-based infinity norm (ADMIN) detector
function [idxhat,bithat] = ADMIN(par,H,y,N0)

% -- preprocessing
% by setting beta to N0/par.Es we get the MMSE estimator in the first iteration
% this is pretty neat as this is a very good detector already
beta = N0/par.Es;%*3; % tweaking this one by 3 improved performance significantly
A = H'*H + beta*eye(par.MT);
L = chol(A,'lower');
yMF = H'*y;

% -- initialization
gamma = (1+sqrt(5))/2;%*2; %% tweaked with 2 to improve performance
alpha = max(real(par.symbols)); % symbol box
zhat = zeros(par.MT,1);
lambda = zeros(par.MT,1);

% -- ADMM loop
for iter=1:par.alg.maxiter
    xhat = (L')\(L\(yMF+beta*(zhat-lambda))); % step 1
    zhat = projinf(par,xhat+lambda,alpha); % step 2
    lambda = lambda-real(gamma*(zhat-xhat)); % step 3
    lambda = real(lambda);
end

% -- hard output detection
[~,idxhat] = min(abs(zhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
bithat = par.bits(idxhat,:);

end




%% Optimized Coordinate Descent (OCD) BOX version
function [idxhat,bithat] = OCDBOX(par,H,y)

% -- initialization
[row, col] = size(H);
alpha = 0; % no regularization for BOX detector
beta = max(real(par.symbols));

% -- preprocessing
dinv = zeros(col,1);
p = zeros(col,1);
for uu=1:col
    normH2 = norm(H(:,uu),2)^2;
    dinv(uu,1) = 1/(normH2+alpha);
    p(uu,1) = dinv(uu)*normH2;
end

r = y;
zold = zeros(col,1);
znew = zeros(col,1);
deltaz = zeros(col,1);

% -- OCD loop
for iters=1:par.alg.maxiter
    for uu=1:col
        tmp = dinv(uu)*(H(:,uu)'*r)+p(uu)*zold(uu);
        znew(uu) = projinf(par,tmp,beta);
        deltaz(uu) = znew(uu)-zold(uu);
        r = r - H(:,uu)*deltaz(uu);
        zold(uu) = znew(uu);
    end
end

[~,idxhat] = min(abs(znew*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
bithat = par.bits(idxhat,:);

end


% project onto alpha infinity-tilde-norm ball
function sproj = projinf(par,s,alpha)
switch par.mod
    case 'BPSK'
        v = real(s);
        sproj = min(abs(v),alpha).*v;
    otherwise
        sr = real(s);
        idxr = abs(sr)>alpha;
        sr(idxr) = sign(sr(idxr))*alpha;
        si = imag(s);
        idxi = abs(si)>alpha;
        si(idxi) = sign(si(idxi))*alpha;
        sproj = sr +1i*si;
end
end

function [idxhat,bithat] = MMSE_approx(par,H,y, sigma2)
n=2;
k=5;

b = H'*y;

omega = 1/(par.MT+par.MR);
inv_P0 = omega*eye(par.MT);
V = sigma2 * inv_P0;
x = inv_P0*b;
R = H'*H;
A = R+sigma2*eye(par.MT);
eigA = eig(A);

omega_exact = 2/(max(eigA) + min(eigA));
for i = 1:n
    x = x + (eye(par.MT)-inv_P0*R-V)^(2^(i-1))*x;
end

for i = (n+1):k
    x = x + omega*(b-R*x);
end
xhat = x;
[~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
bithat = par.bits(idxhat,:);
end

% function [idxhat,bithat] = bMMSE(par,H,y,N0)
% xhat = (H'*H+(N0/par.Es)*eye(par.MT))\(H'*y);
% [~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
% bithat = par.bits(idxhat,:);
% end
function [idxhat,bithat] = MMSE_approx2(par,H,y, N0)
n=1;
k=5;

omega = 1/(par.MT+par.MR);
A = H'*H +(N0/par.Es)*eye(par.MT);
b = H'*y;
eta = par.MR/par.MT;
alpha = 2*eta*(1+eta)/(par.MR*(1+eta^2));
phi = -eta^2/(par.MR^2*(1+eta^2));
X = alpha*eye(par.MT,par.MT) + phi*A;

eigA = eig(A);

omega_exact = 2/(max(eigA) + min(eigA));
% for i = 1 :2
%    X = X*(2*eye(par.MT,par.MT) - A*X); 
% end
% s = X*b;


s0 = X*b;
s1 = 2*s0 - X*A*s0;

s = s1;
for i = 1:3
    s = s + omega*(b-A*s);
end


xhat = s;
[~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
bithat = par.bits(idxhat,:);
end

function [idxhat,bithat] = MMSE_approx4(par,H,y, N0)
n=1;
k=5;

omega = 1/(par.MT+par.MR);
A = H'*H +(N0/par.Es)*eye(par.MT);
b = H'*y;
D = diag(diag(A));
Dinv = inv(D);
eta = par.MR/par.MT;
alpha = 2*eta*(1+eta)/(par.MR*(1+eta^2));
phi = -eta^2/(par.MR^2*(1+eta^2));
X = alpha*eye(par.MT,par.MT) + phi*A;

% for i = 1 :2
%    X = X*(2*eye(par.MT,par.MT) - A*X); 
% end
% s = X*b;
stair_matrix = stair(A);
INV_stair = inv(stair_matrix);

s0 = X*b;
s1 = 2*s0 - X*A*s0;

s = s1;
for i = 1:2
    s =  INV_stair*((stair_matrix-A)*s +b);
end


xhat = s;
[~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
bithat = par.bits(idxhat,:);
end
function [idxhat,bithat] = MMSE_approx3(par,H,y, N0)

A = H'*H +(N0/par.Es)*eye(par.MT);
b = H'*y;
D = diag(diag(A));
Dinv = inv(D);
x = Dinv*b;
r = b - A*x;
stair_matrix = stair(A);
INV_stair = inv(stair_matrix);
% eta = par.MR/par.MT;
% alpha = 2*eta*(1+eta)/(par.MR*(1+eta^2));
% phi = -eta^2/(par.MR^2*(1+eta^2));
% X = alpha*eye(par.MT,par.MT) + phi*A;

% for i = 1 :2
%    X = X*(2*eye(par.MT,par.MT) - A*X); 
% end
% s = X*b;


% s0 = X*b;
% s1 = 2*s0 - X*A*s0;
% s = s1;
p = A*r;
u = r'*r/(p'*r);
s = x + u*r + Dinv*(r - u*p);
for i = 1:2
    s =  INV_stair*((stair_matrix-A)*s +b);
end


xhat = s;
[~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
bithat = par.bits(idxhat,:);
end

function [idxhat,bithat] = MMSE_approx_exact_eigen(par,H,y, N0)
n=1;
k=5;

omega = 1/(par.MT+par.MR);
A = H'*H +(N0/par.Es)*eye(par.MT);
eigA = eig(A);
omega_exact = 2/(max(eigA) + min(eigA));
lamda_mid = 1/omega_exact;
b = H'*y;
eta = par.MR/par.MT;

alpha = 4*lamda_mid/(2*lamda_mid^2 -(max(eigA) -lamda_mid)^2 );
phi = -2/(2*lamda_mid^2 -(max(eigA) -lamda_mid)^2 );
% alpha = 2*eta*(1+eta)/(par.MR*(1+eta^2));
% phi = -eta^2/(par.MR^2*(1+eta^2));
X = alpha*eye(par.MT,par.MT) + phi*A;


% for i = 1 :2
%    X = X*(2*eye(par.MT,par.MT) - A*X); 
% end
% s = X*b;


s0 = X*b;
s1 = 2*s0 - X*A*s0;

s = s1;
for i = 1:3
    s = s + omega_exact*(b-A*s);
end


xhat = s;
[~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
bithat = par.bits(idxhat,:);
end
function [idxhat,bithat] = MMSE_approx2_real(par,Hc,yc, N0)
H = [real(Hc) -imag(Hc);
    imag(Hc)  real(Hc)];
y = [real(yc)' imag(yc)']';
n=1;
k=5;

omega = 1/(par.MT+par.MR);
A = H'*H +(N0/par.Es)*eye(2*par.MT);
b = H'*y;
eta = par.MR/par.MT;
alpha = 2*eta*(1+eta)/(par.MR*(1+eta^2));
phi = -eta^2/(par.MR^2*(1+eta^2));
X = alpha*eye(2*par.MT) + phi*A;

eigA = eig(A);
% for i = 1 :2
%    X = X*(2*eye(par.MT,par.MT) - A*X); 
% end
% s = X*b;


s0 = X*b;
s1 = 2*s0 - X*A*s0;

s = s1;
for i = 1:3
    s = s + omega*(b-A*s);
end


xhat = s;
% [~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
% bithat = par.bits(idxhat,:);

[~,idxhat] = min(abs(xhat*ones(1,length(par.RealSymbols))-ones(2*par.MT,1)*par.RealSymbols).^2,[],2);
% x = par.RealSymbols(idxhat)';
idxhat1 = zeros(2*par.MT,1);
for i = 1:par.MT
    idxhat1(2*i-1,1) = idxhat(i,1);
    idxhat1(2*i,1)   = idxhat(i+par.MT,1);
end


bithat = par.realbits(idxhat1,:);
bithat = reshape(bithat', par.Q, par.MT)';
end
