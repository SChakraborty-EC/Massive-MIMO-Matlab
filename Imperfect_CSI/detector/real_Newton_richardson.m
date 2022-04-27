function xhat = real_Newton_richardson(par,H,y, N0)

A = H'*H +(N0/par.Es)*eye(2*par.MT);
b = H'*y;
eta = par.MR/par.MT;
alpha = 2*eta*(1+eta)/(par.MR*(1+eta^2));
phi = -eta^2/(par.MR^2*(1+eta^2));
omega = 1/(par.MT+par.MR);
X = alpha*eye(2*par.MT,2*par.MT) + phi*A;

s0 = X*b;
s1 = 2*s0 - X*A*s0;

s = s1;
for i = 1:par.alg.maxiter-1
    s = s + omega*(b-A*s);
end

xhat = s;
% [~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
% bithat = par.bits(idxhat,:);
end