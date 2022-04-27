function x_hat = SOR(par,H,y,N0)
A = H'*H+(N0/par.Es)*eye(par.MT);

D = diag(diag(A));

U = triu(A,1);
L = tril(A,-1);


yMF = H'*y;

inv_D  = inv(D);
x_hat =  inv(D)*yMF;
%  x_hat = inv_D*yMF;
% x_hat = zeros(par.MT,1);
omega = 1.05; 
inv_omega = 1/omega;

for k = 1:par.alg.maxiter
x_hat = inv(D*inv_omega +L) *(yMF+((inv_omega-1)*D-U)*x_hat);

%     xhat = inv(D+L)*(yMF - U*xhat);
end


% [~,idxhat] = min(abs(x_hat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
% bithat = par.bits(idxhat,:);