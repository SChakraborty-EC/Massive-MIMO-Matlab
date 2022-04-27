%% Conjugate Gradient massive MIMO detection

function bithat = CG_imperfect_CSI(par,H,yBlock,N0)
    A = H'*H+(N0/par.Es)*eye(par.MT);
bithat = zeros(par.MT,par.Q,par.block_size);
for i=1:par.block_size
    y = yBlock(:,i);

    MF = H'*y;
    
    r = MF;
    p = r;
    v = zeros(par.MT,1);
    
    for k = 1:par.alg.maxiter
        e = A*p;
        alpha  = norm(r)^2/(p'*e);
        v = v+alpha*p;
        new_r = r-alpha*e;
        beta = norm(new_r)^2/norm(r)^2;
        p = new_r+beta*p;
        r = new_r;
    end
    
    xhat = v;
    [~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
    bithat(:,:,i) =  par.bits(idxhat,:);   
end