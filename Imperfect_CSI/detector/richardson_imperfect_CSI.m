function bithat = richardson_imperfect_CSI(par,H,yBlock,N0)
A = H'*H+(N0/par.Es)*eye(par.MT);


bithat = zeros(par.MT,par.Q,par.block_size);
%detect the symbols using inperfect CSI
for i=1:par.block_size
    y = yBlock(:,i);
    
    yMF = H'*y;
    D = diag(diag(A));
    inv_D  = inv(D);
    omega = 1/(par.MT + par.MR);
    xhat = inv_D*yMF;
    
    for k = 1:par.alg.maxiter
        xhat = xhat + omega*(yMF - A*xhat);
    end
    [~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
    bithat(:,:,i) =  par.bits(idxhat,:);
end

% [~,idxhat] = min(abs(x_hat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
% bithat = par.bits(idxhat,:);