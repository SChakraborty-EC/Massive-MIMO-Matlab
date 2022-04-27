function xhat = real_MMSE_approx(par,H,y, sigma2)
n=2;


b = H'*y;

omega = 1/(par.MT+par.MR);
inv_P0 = omega*eye(2*par.MT);
V = sigma2 * inv_P0;
x = inv_P0*b;
R = H'*H;
% A = R+sigma2*eye(par.MT);

for i = 1:n
    x = x + (eye(2*par.MT)-inv_P0*R-V)^(2^(i-1))*x;
end

for i = (n+1):par.alg.maxiter
    x = x + omega*(b-R*x);
end
xhat = x;
% [~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
% bithat = par.bits(idxhat,:);
end