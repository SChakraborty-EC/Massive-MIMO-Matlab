function x_hat = jacobi(par,H,y,N0)
A = H'*H+(N0/par.Es)*eye(par.MT);

yMF = H'*y;
D = diag(diag(A));
inv_D  = inv(D);

x_hat = inv_D*yMF;

for k = 1:par.alg.maxiter
x_hat = inv_D*(yMF + (D-A)*x_hat);   
end


% [~,idxhat] = min(abs(x_hat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
% bithat = par.bits(idxhat,:);