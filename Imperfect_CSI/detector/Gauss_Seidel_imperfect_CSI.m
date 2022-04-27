%% Gauss-Seidel massive MIMO detection

function bithat = Gauss_Seidel_imperfect_CSI(par,H,yBlock,N0)
A = H'*H+(N0/par.Es)*eye(par.MT);
D = diag(diag(A));
U = triu(A,1);
L = tril(A,-1);
bithat = zeros(par.MT,par.Q,par.block_size);

%detect the symbols using inperfect CSI
for i=1:par.block_size
    y = yBlock(:,i);
    
    MF = H'*y;
    % xhat = diag(inv(D));% inv(D)*MF;  %%% Check Gauss Seidel detection paper
    xhat =  inv(D)*MF;
    for j = 1:par.alg.maxiter
        xhat = inv(D+L)*(MF - U*xhat);
    end
    [~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
    bithat(:,:,i) =  par.bits(idxhat,:);
end