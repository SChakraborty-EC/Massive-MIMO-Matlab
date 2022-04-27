function bithat = chebyshev_imperfect_CSI(par,H,yBlock,N0)
lambda_max = par.MR*(1+sqrt(par.MT/par.MR))^2;
lambda_min = par.MR*(1-sqrt(par.MT/par.MR))^2;

alpha = (lambda_max + lambda_min)/(lambda_max-lambda_min);
beta  = (lambda_max + lambda_min)/2;

bithat = zeros(par.MT,par.Q,par.block_size);
%detect the symbols using inperfect CSI
for i=1:par.block_size
y = yBlock(:,i);




yMF = H'*y;
roh = 2*alpha^2/((2*alpha^2-1)*beta);
phi = 1/(2*alpha^2 -1);
s_hat = yMF/beta;
h = H*s_hat;
r = yMF - (N0/par.Es)*s_hat - H'*h;
sigma = r/beta;

for k = 1:par.alg.maxiter
   s_hat = s_hat + sigma;
   h = H*s_hat;
   r = yMF - (N0/par.Es)*s_hat - H'*h;
   sigma = roh*r + phi*sigma;
   roh_old = roh;
   roh = 4*alpha^2/(4*(alpha^2)*beta -(beta^2)*roh_old);
   phi = beta^2*roh_old/(4*(alpha^2)*beta -(beta^2)*roh_old);    
end

[~,idxhat] = min(abs(s_hat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
bithat(:,:,i) =  par.bits(idxhat,:);
end
% [~,idxhat] = min(abs(s*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
% bithat = par.bits(idxhat,:);