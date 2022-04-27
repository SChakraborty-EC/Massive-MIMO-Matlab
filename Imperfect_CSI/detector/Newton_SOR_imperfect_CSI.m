function bithat = Newton_SOR_imperfect_CSI(par,H,yBlock,N0)
A = H'*H+(N0/par.Es)*eye(par.MT);

D = diag(diag(A));

U = triu(A,1);
L = tril(A,-1);
inv_D  = inv(D);
bithat = zeros(par.MT,par.Q,par.block_size);
%detect the symbols using inperfect CSI
for i=1:par.block_size
    y = yBlock(:,i);
    
    yMF = H'*y;
       
    eta = par.MR/par.MT;
    alpha = 2*eta*(1+eta)/(par.MR*(1+eta^2));
    phi = -eta^2/(par.MR^2*(1+eta^2));
    omega = 1/(par.MT+par.MR);
    X = alpha*eye(par.MT,par.MT) + phi*A;
    
    s0 = X*yMF;
    xhat = 2*s0 - X*A*s0;

    % x_hat =  inv(D)*yMF;
    %  x_hat = inv_D*yMF;
    % x_hat = zeros(par.MT,1);
    %omega = 1.05; % for SOR
    omega =1;  % for GS
    inv_omega = 1/omega;
    
    for k = 1:par.alg.maxiter-1
        xhat = inv(D*inv_omega +L) *(yMF+((inv_omega-1)*D-U)*xhat);
        
        %     xhat = inv(D+L)*(yMF - U*xhat);
    end
    [~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
    bithat(:,:,i) =  par.bits(idxhat,:);
end