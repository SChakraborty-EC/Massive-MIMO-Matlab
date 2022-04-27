function xhat = Newton_chebyshev(par,H,y, N0)
% --- parameters for newton iteration
A = H'*H +(N0/par.Es)*eye(par.MT);
b = H'*y;
eta = par.MR/par.MT;
alpha = 2*eta*(1+eta)/(par.MR*(1+eta^2));
phi = -eta^2/(par.MR^2*(1+eta^2));
X = alpha*eye(par.MT) + phi*A;

s0 = X*b;
s1 = 2*s0 - X*A*s0;
s = s1;


% ----- Parameters for chebyshev
lambda_max = par.MR*(1+sqrt(par.MT/par.MR))^2;
lambda_min = par.MR*(1-sqrt(par.MT/par.MR))^2;

alpha = (lambda_max+lambda_min)/(lambda_max - lambda_min);
beta  = (lambda_max + lambda_min)/2;

yMF = H'*y;

roh = 2*alpha^2/((2*alpha^2-1)*beta);
phi = 1/(2*alpha^2 -1);

% s = yMF/beta;
h = H*s;
r = yMF - (N0/par.Es)*s - H'*h;
sigma = r/beta;

for k = 1:par.alg.maxiter-1
s = s + sigma;
h = H*s;
r = yMF -(N0/par.Es)*s - H'*h;
sigma = roh * r + phi*sigma;
phi = (beta^2)*roh/(4*alpha^2*beta - beta^2*roh);
roh = 4*alpha^2/(4*alpha^2*beta - (beta^2)*roh);
end

xhat = s;

% [~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
% bithat = par.bits(idxhat,:);