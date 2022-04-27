%% Optimized Coordinate Descent (OCD) BOX version
function znew = OCDBOX(par,H,y)

% -- initialization
[row, col] = size(H);
alpha = 0; % no regularization for BOX detector
beta = max(real(par.symbols));

% -- preprocessing
dinv = zeros(col,1);
p = zeros(col,1);
for uu=1:col
    normH2 = norm(H(:,uu),2)^2;
    dinv(uu,1) = 1/(normH2+alpha);
    p(uu,1) = dinv(uu)*normH2;
end

r = y;
zold = zeros(col,1);
znew = zeros(col,1);
deltaz = zeros(col,1);

% -- OCD loop
for iters=1:par.alg.maxiter
    for uu=1:col
        tmp = dinv(uu)*(H(:,uu)'*r)+p(uu)*zold(uu);
        znew(uu) = projinf(par,tmp,beta);
        deltaz(uu) = znew(uu)-zold(uu);
        r = r - H(:,uu)*deltaz(uu);
        zold(uu) = znew(uu);
    end
end

% [~,idxhat] = min(abs(znew*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
% bithat = par.bits(idxhat,:);

end

% project onto alpha infinity-tilde-norm ball
function sproj = projinf(par,s,alpha)
switch par.mod
    case 'BPSK'
        v = real(s);
        sproj = min(abs(v),alpha).*v;
    otherwise
        sr = real(s);
        idxr = abs(sr)>alpha;
        sr(idxr) = sign(sr(idxr))*alpha;
        si = imag(s);
        idxi = abs(si)>alpha;
        si(idxi) = sign(si(idxi))*alpha;
        sproj = sr +1i*si;
end
end