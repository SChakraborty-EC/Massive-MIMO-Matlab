function x = MMSE_init(par,H,y,N0)
switch (par.MMSE_init) % select algorithms
    case 'ZF', % zero-forcing detection
        xhat = ZF(par,H,y);        
    case 'bMMSE', % biased MMSE detector
        xhat = real_bMMSE(par,H,y,N0);        
    case 'ML', % ML detection using sphere decoding
        xhat = ML(par,H,y);
    case 'SIMO' % SIMO lower bound
        xhat = SIMO(par,H,y,N0,s);               
    case 'Gauss-Seidel' % Gauss-Seidel detector
        xhat = real_Gauss_Seidel(par,H,y,N0);        
    case 'ADMIN' % ADMM-based Infinity Norm detector
        xhat = real_ADMIN(par,H,y,N0);        
    case 'OCDBOX' % co-ordinate descent (optimized) detector
        xhat = real_OCDBOX(par,H,y);        
    case 'Neumann' % coordinate descent
        xhat = real_Neumann(par,H,y,N0);        
    case 'Conjugate-Gradient' % conjugate gradient detector
        xhat = real_CG(par,H,y,N0);        
    case 'richardson' % conjugate gradient detector
        xhat = real_richardson(par,H,y,N0);        
    case 'SOR' % conjugate gradient detector
        xhat = real_SOR(par,H,y,N0);        
    case 'Newton_iteration' % conjugate gradient detector
        xhat = real_Newton_iteration(par,H,y,N0);        
    case 'jacobi' % conjugate gradient detector
        xhat = real_jacobi(par,H,y,N0);  
    case 'chebyshev' % conjugate gradient detector
        xhat = real_chebyshev(par,H,y,N0); 
    case 'Newton_richardson' % conjugate gradient detector
        xhat = real_Newton_richardson(par,H,y,N0);       
    case 'Newton_chebyshev' % conjugate gradient detector
        xhat = real_Newton_chebyshev(par,H,y,N0);     
    case 'Newton_GS' % conjugate gradient detector
        xhat = real_Newton_GS(par,H,y,N0);  
    case 'MMSEapprox', % biased MMSE detector
        xhat = real_MMSE_approx(par,H,y, N0);        
    otherwise,
        error('par.detector type not defined.')
end


[~,idxhat] = min(abs(xhat*ones(1,length(par.RealSymbols))-ones(2*par.MT,1)*par.RealSymbols).^2,[],2);
x = par.RealSymbols(idxhat)';

% [~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
% x = par.symbols(idxhat)';
end


% function [idxhat,bithat] = MF(par,H,y)
%
% xhat = H' * y / norm(H(:));
% [~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
% bithat = par.bits(idxhat,:);
% end

% %% MMSE detector (MMSE)
% function xhat = MMSE(par,H,y,N0)
% xhat = (H'*H+(N0/par.Es)*eye(2*par.MT))\(H'*y);
% % [~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
% % bithat = par.bits(idxhat,:);
% end
%
% %% SIMO bound
% function [idxhat,bithat] = SIMO(par,H,y,s)
% z = y-H*s;
% xhat = zeros(par.MT,1);
% for m=1:par.MT
%     hm = H(:,m);
%     yhat = z+hm*s(m,1);
%     xhat(m,1) = hm'*yhat/norm(hm,2)^2;
% end
% [~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
% bithat = par.bits(idxhat,:);
% end
%
%
% %% Neumann-Series Approximation based massive MIMO detection
%
% function xhat = Neumann(par,H,y,N0)
% A = H'*H+(N0/par.Es)*eye(2*par.MT);
% MF = H'*y;
%
% D = diag(diag(A));
% E = triu(A,1)+tril(A,-1);
% Ainv = 0;
% for i = 0:par.alg.maxiter
%     Ainv = Ainv+((-inv(D)*E)^i)*inv(D);
% end
%
% xhat = Ainv*MF;
% % [~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
% % bithat = par.bits(idxhat,:);
% end
%
% %% Gauss-Seidel massive MIMO detection
%
% function xhat = Gauss_Seidel(par,H,y,N0)
% A = H'*H+(N0/par.Es)*eye(2*par.MT);
% MF = H'*y;
%
% D = diag(diag(A));
% E = -triu(A,1);
% F = -tril(A,-1);
%
% xhat = diag(inv(D));% inv(D)*MF;  %%% Check Gauss Seidel detection paper
% for i = 0:par.alg.maxiter
%     xhat = inv(D-E)*(F*xhat+MF);
% end
%
% % [~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
% % bithat = par.bits(idxhat,:);
% end
%
%
%
% %% Conjugate Gradient massive MIMO detection
%
% function xhat = CG(par,H,y,N0)
% A = H'*H+(N0/par.Es)*eye(2*par.MT);
% MF = H'*y;
%
% r = MF;
% p = r;
% v = zeros(2*par.MT,1);
%
% for k = 1:par.alg.maxiter
%     e = A*p;
%     alpha  = norm(r)^2/(p'*e);
%     v = v+alpha*p;
%     new_r = r-alpha*e;
%     beta = norm(new_r)^2/norm(r)^2;
%     p = new_r+beta*p;
%     r = new_r;
% end
%
% xhat = v;
%
% % [~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
% % bithat = par.bits(idxhat,:);
%
% end
%
%
% %% ADMM-based infinity norm (ADMIN) detector
% function [idxhat,bithat] = ADMIN(par,H,y,N0)
%
% % -- preprocessing
% % by setting beta to N0/par.Es we get the MMSE estimator in the first iteration
% % this is pretty neat as this is a very good detector already
% beta = N0/par.Es;%*3; % tweaking this one by 3 improved performance significantly
% A = H'*H + beta*eye(2*par.MT);
% L = chol(A,'lower');
% yMF = H'*y;
%
% % -- initialization
% gamma = (1+sqrt(5))/2;%*2; %% tweaked with 2 to improve performance
% alpha = max(real(par.symbols)); % symbol box
% zhat = zeros(par.MT,1);
% lambda = zeros(par.MT,1);
%
% % -- ADMM loop
% for iter=1:par.alg.maxiter
%     xhat = (L')\(L\(yMF+beta*(zhat-lambda))); % step 1
%     zhat = projinf(par,xhat+lambda,alpha); % step 2
%     lambda = lambda-real(gamma*(zhat-xhat)); % step 3
%     lambda = real(lambda);
% end
%
% % -- hard output detection
% [~,idxhat] = min(abs(zhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
% bithat = par.bits(idxhat,:);
%
% end
%
%
%
%
% %% Optimized Coordinate Descent (OCD) BOX version
% function [idxhat,bithat] = OCDBOX(par,H,y)
%
% % -- initialization
% [row, col] = size(H);
% alpha = 0; % no regularization for BOX detector
% beta = max(real(par.symbols));
%
% % -- preprocessing
% dinv = zeros(col,1);
% p = zeros(col,1);
% for uu=1:col
%     normH2 = norm(H(:,uu),2)^2;
%     dinv(uu,1) = 1/(normH2+alpha);
%     p(uu,1) = dinv(uu)*normH2;
% end
%
% r = y;
% zold = zeros(col,1);
% znew = zeros(col,1);
% deltaz = zeros(col,1);
%
% % -- OCD loop
% for iters=1:par.alg.maxiter
%     for uu=1:col
%         tmp = dinv(uu)*(H(:,uu)'*r)+p(uu)*zold(uu);
%         znew(uu) = projinf(par,tmp,beta);
%         deltaz(uu) = znew(uu)-zold(uu);
%         r = r - H(:,uu)*deltaz(uu);
%         zold(uu) = znew(uu);
%     end
% end
%
% [~,idxhat] = min(abs(znew*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
% bithat = par.bits(idxhat,:);
%
% end
%
%
% % project onto alpha infinity-tilde-norm ball
% function sproj = projinf(par,s,alpha)
% switch par.mod
%     case 'BPSK'
%         v = real(s);
%         sproj = min(abs(v),alpha).*v;
%     otherwise
%         sr = real(s);
%         idxr = abs(sr)>alpha;
%         sr(idxr) = sign(sr(idxr))*alpha;
%         si = imag(s);
%         idxi = abs(si)>alpha;
%         si(idxi) = sign(si(idxi))*alpha;
%         sproj = sr +1i*si;
% end
% end
%
% function xhat = MMSE_approx(par,H,y, sigma2)
% n=2;
% k=4;
%
% b = H'*y;
%
% omega = 1/(2*par.MT+2*par.MR);
% inv_P0 = omega*eye(2*par.MT);
% V = sigma2 * inv_P0;
% x = inv_P0*b;
% R = H'*H;
%
% for i = 1:n
%     x = x + (eye(2*par.MT)-inv_P0*R-V)^(2^(i-1))*x;
% end
%
% for i = (n+1):k
%     x = x + omega*(b-R*x);
% end
% xhat = x;
% % [~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
% % bithat = par.bits(idxhat,:);
% end
%
% function xhat = bMMSE(par,H,y,N0)
% xhat = (H'*H+(N0/par.Es)*eye(2*par.MT))\(H'*y);
% % [~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
% % bithat = par.bits(idxhat,:);
% end
%
% function xhat = MMSE_approx2_real(par,H,y, N0)
% % H = [real(Hc) -imag(Hc);
% %     imag(Hc)  real(Hc)];
% % y = [real(yc)' imag(yc)']';
% n=1;
% k=5;
%
% omega = 1/(par.MT+par.MR);
% A = H'*H +(N0/par.Es)*eye(2*par.MT);
% b = H'*y;
% eta = par.MR/par.MT;
% alpha = 2*eta*(1+eta)/(par.MR*(1+eta^2));
% phi = -eta^2/(par.MR^2*(1+eta^2));
% X = alpha*eye(2*par.MT) + phi*A;
%
% eigA = eig(A);
% % for i = 1 :2
% %    X = X*(2*eye(par.MT,par.MT) - A*X);
% % end
% % s = X*b;
%
%
% s0 = X*b;
% s1 = 2*s0 - X*A*s0;
%
% s = s1;
% for i = 1:3
%     s = s + omega*(b-A*s);
% end
%
%
% xhat = s;
% % [~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
% % bithat = par.bits(idxhat,:);
%
% % [~,idxhat] = min(abs(xhat*ones(1,length(par.RealSymbols))-ones(2*par.MT,1)*par.RealSymbols).^2,[],2);
% % % x = par.RealSymbols(idxhat)';
% % idxhat1 = zeros(2*par.MT,1);
% % for i = 1:par.MT
% %     idxhat1(2*i-1,1) = idxhat(i,1);
% %     idxhat1(2*i,1)   = idxhat(i+par.MT,1);
% % end
% %
% %
% % bithat = par.realbits(idxhat1,:);
% % bithat = reshape(bithat', par.Q, par.MT)';
% end
%
