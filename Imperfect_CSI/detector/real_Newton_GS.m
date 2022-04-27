function xhat = real_Newton_GS(par,H,y, N0)

A = H'*H +(N0/par.Es)*eye(2*par.MT);
b = H'*y;
eta = par.MR/par.MT;
alpha = 2*eta*(1+eta)/(par.MR*(1+eta^2));
phi = -eta^2/(par.MR^2*(1+eta^2));
X = alpha*eye(2*par.MT,2*par.MT) + phi*A;

s0 = X*b;
xhat = 2*s0 - X*A*s0;


D = diag(diag(A));
E = -triu(A,1);
F = -tril(A,-1);

for k = 1:par.alg.maxiter-1
    xhat = inv(D-E)*(F*xhat+b);
end


% [~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
% bithat = par.bits(idxhat,:);
end