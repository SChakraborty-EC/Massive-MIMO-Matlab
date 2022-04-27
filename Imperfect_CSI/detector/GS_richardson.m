function xhat = GS_richardson(par,H,y,N0)
A = H'*H+(N0/par.Es)*eye(par.MT);
MF = H'*y;
omega = 1/(par.MT+par.MR);

D = diag(diag(A));
U = triu(A,1);
L = tril(A,-1);

% xhat = diag(inv(D));% inv(D)*MF;  %%% Check Gauss Seidel detection paper
xhat =  inv(D)*MF;
for i = 1:1
    xhat = inv(D+L)*(MF - U*xhat);
end

for i = 1:par.alg.maxiter-1
    xhat = xhat + omega*(MF-A*xhat);
end
% [~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
% bithat = par.bits(idxhat,:);
end