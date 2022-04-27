%% biased MMSE detector (bMMSE)
function xhat = bMMSE(par,H,y,N0)
xhat = (H'*H+(N0/par.Es)*eye(par.MT))\(H'*y);
% [~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
% bithat = par.bits(idxhat,:);
end