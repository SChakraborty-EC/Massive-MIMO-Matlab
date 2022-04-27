%% Neumann-Series Approximation based massive MIMO detection

function xhat = real_Neumann(par,H,y,N0)
A = H'*H+(N0/par.Es)*eye(2*par.MT);
MF = H'*y;

D = diag(diag(A));
E = triu(A,1)+tril(A,-1);
Ainv = 0;
for i = 1:par.alg.maxiter
    Ainv = Ainv+((-inv(D)*E)^i)*inv(D);
end

xhat = Ainv*MF;
% [~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
% bithat = par.bits(idxhat,:);
end