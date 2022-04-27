%% biased MMSE detector (bMMSE)
function xhat = real_bMMSE(par,H,y,N0)
xhat = (H'*H+(N0/par.Es)*eye(2*par.MT))\(H'*y);
% [~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
% bithat = par.bits(idxhat,:);
end