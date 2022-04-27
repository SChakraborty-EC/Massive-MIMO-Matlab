function bithat = Newton_GS_imperfect_CSI(par,H,yBlock, N0)

A = H'*H +(N0/par.Es)*eye(par.MT);
D = diag(diag(A));
E = -triu(A,1);
F = -tril(A,-1);

eta = par.MR/par.MT;
alpha = 2*eta*(1+eta)/(par.MR*(1+eta^2));
phi = -eta^2/(par.MR^2*(1+eta^2));
X = alpha*eye(par.MT,par.MT) + phi*A;

bithat = zeros(par.MT,par.Q,par.block_size);
%detect the symbols using inperfect CSI
for i=1:par.block_size
y = yBlock(:,i);

b = H'*y;


s0 = X*b;
xhat = 2*s0 - X*A*s0;

for k = 1:par.alg.maxiter-1
    xhat = inv(D-E)*(F*xhat+b);
end

[~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
bithat(:,:,i) =  par.bits(idxhat,:);
end