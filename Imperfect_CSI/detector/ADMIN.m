%% ADMM-based infinity norm (ADMIN) detector
function zhat = ADMIN(par,H,y,N0)

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
% [~,idxhat] = min(abs(zhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
% bithat = par.bits(idxhat,:);

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