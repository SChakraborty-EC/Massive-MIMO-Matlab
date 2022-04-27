%% Gauss-Seidel massive MIMO detection

function xhat = real_Gauss_Seidel(par,H,y,N0)
A = H'*H+(N0/par.Es)*eye(2*par.MT);
MF = H'*y;

D = diag(diag(A));
E = -triu(A,1);
F = -tril(A,-1);

xhat = diag(inv(D));% inv(D)*MF;  %%% Check Gauss Seidel detection paper
for i = 1:par.alg.maxiter
    xhat = inv(D-E)*(F*xhat+MF);
end

% [~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
% bithat = par.bits(idxhat,:);
end