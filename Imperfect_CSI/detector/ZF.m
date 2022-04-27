%% zero-forcing (ZF) detector
function x_hat = ZF(par,H,y)
xhat = H\y;
% [~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
% bithat = par.bits(idxhat,:);
end