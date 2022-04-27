%% Matched filter

function xhat = MF(par,H,y)

xhat = H' * y / norm(H(:));
% [~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
% bithat = par.bits(idxhat,:);
end