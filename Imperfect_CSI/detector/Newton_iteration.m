function xhat = Newton_iteration(par,H,y, N0)

A = H'*H +(N0/par.Es)*eye(par.MT);
b = H'*y;
eta = par.MR/par.MT;
alpha = 2*eta*(1+eta)/(par.MR*(1+eta^2));
phi = -eta^2/(par.MR^2*(1+eta^2));

X = alpha*eye(par.MT,par.MT) + phi*A;

for i = 1 :2
   X = X*(2*eye(par.MT,par.MT) - A*X);
end
s = X*b;

xhat = s;
% [~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
% bithat = par.bits(idxhat,:);
end