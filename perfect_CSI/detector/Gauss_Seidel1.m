%% Gauss-Seidel massive MIMO detection

function xhat = Gauss_Seidel1(par,H,y,N0)
A = H'*H+(N0/par.Es)*eye(par.MT);
MF = H'*y;

D = diag(diag(A));
U = triu(A,1);
L = tril(A,-1);

% xhat = diag(inv(D));% inv(D)*MF;  %%% Check Gauss Seidel detection paper
xhat =  inv(D)*MF;
for i = 1:par.alg.maxiter
    xhat = inv(D+L)*(MF - U*xhat);
end

% [~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
% bithat = par.bits(idxhat,:);
end