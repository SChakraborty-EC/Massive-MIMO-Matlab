function x_hat = steepest_jacobi(par,H,y,N0)
G = H'*H;
A = G+(N0/par.Es)*eye(par.MT);

b = H'*y;
D = diag(diag(A));
inv_D  = inv(D);
x_hat = inv_D*b;
r0 = b - A*x_hat;

%% first iteration

p = A*r0;
u = (r0'*r0)/(p'*r0);
x_hat = x_hat + u*r0 + inv_D*(r0 - u*p);

%% jacobi iteration
for k = 1:par.alg.maxiter-1
x_hat = inv_D*((D-A)*x_hat + b);
end



% [~,idxhat] = min(abs(x_hat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
% bithat = par.bits(idxhat,:);