%% SIMO bound
function bithat = SIMO(par,H,y,s)
z = y-H*s;
xhat = zeros(par.MT,1);
for m=1:par.MT
    hm = H(:,m);
    yhat = z+hm*s(m,1);
    xhat(m,1) = hm'*yhat/norm(hm,2)^2;
end
[~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
bithat = par.bits(idxhat,:);
end