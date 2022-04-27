%% MMSE detector (MMSE)
function bithat = MMSE_imperfect_CSI(par,H,yBlock,N0)


bithat = zeros(par.MT,par.Q,par.block_size);
%detect the symbols using inperfect CSI
for i=1:par.block_size
y = yBlock(:,i);

xhat = (H'*H+(N0/par.Es)*eye(par.MT))\(H'*y);
[~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
bithat(:,:,i) =  par.bits(idxhat,:);
end