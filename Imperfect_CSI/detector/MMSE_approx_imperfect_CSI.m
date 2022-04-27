function bithat = MMSE_approx_imperfect_CSI(par,H,yBlock, sigma2)
n=2;

bithat = zeros(par.MT,par.Q,par.block_size);
%detect the symbols using inperfect CSI
for i=1:par.block_size
    y = yBlock(:,i);
    
    b = H'*y;
    temp = eye(par.MT, par.MT);
    
    omega = 1/(par.MT+par.MR);
    inv_P0 = omega*temp;
    V = sigma2 * inv_P0;
    x = inv_P0*b;
    R = H'*H;
    % A = R+sigma2*eye(par.MT);
    
    x = x + (temp-inv_P0*(H'*H)-V)*x;
    x = x + ((temp-inv_P0*(H'*H)-V)^2) *x;
    
    for j = (n+1):par.alg.maxiter
        x = x + omega*(b-R*x);
    end
    xhat = x;
    [~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
    bithat(:,:,i) =  par.bits(idxhat,:);
end
end