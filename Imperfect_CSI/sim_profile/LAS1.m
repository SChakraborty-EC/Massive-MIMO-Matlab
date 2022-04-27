function bithat = LAS1(par, Hc,yc, N0)
H = zeros(2*par.MR,2*par.MT);
%% Real valued decomposition for channel matrix
for i = 1:par.MR % ROW
    for j = 1:par.MT % column
        H(2*i-1, 2*j-1) =  real(Hc(i,j));
        H(2*i-1, 2*j)   = -imag(Hc(i,j));
        H(2*i, 2*j-1)   =  imag(Hc(i,j));
        H(2*i, 2*j)     =  real(Hc(i,j));
    end
end
y = zeros(2*par.MR,1);
for i = 1:par.MR % ROW
    y(2*i-1,1) = real(yc(i,1));
    y(2*i,1)   = imag(yc(i,1));
end

% Get the initial solution
k = 1;
% x(:,r) = ZF_init(par, H,y);
d(:,k) = MMSE_init(par,H,y,N0);
% d = zeros(2*par.MT,1);
% for i = 1:par.MT % ROW
%     d(2*i-1,k) = real(s_hat(i,1));
%     d(2*i,k)   = imag(s_hat(i,1));
% end
G = H'*H;
% G=bandmatrix(G,par.w); 
if par.band_approx == 1
   G=bandmatrix_mod(G,par.w);
end

z = H'*y - H'*H* d(:,k);

terminate =0;

lopt = zeros(1, 2*par.MT);
tilde_lopt = zeros(1, 2*par.MT);
d_tilde = zeros( 2*par.MT, 1);
cost = zeros(1, 2*par.MT);
k = 1;

while terminate ~= 1
    
    for i = 1:2*par.MT
        
        lopt(1,i) = 2* round(abs(z(i,1))/(2*G(i,i)));
        d_tilde (i) = d(i,k) + lopt(1,i)* sign(z(i));
        
        if d_tilde (i,1) > par.M
            d_tilde (i,1) = par.M;
            
        elseif d_tilde (i,1) < -par.M
            d_tilde (i,1) = -par.M;
        else
            d_tilde (i,1) = d_tilde (i,1);
        end
        tilde_lopt(1,i) = d_tilde (i,1) - d(i,k);
        cost(1,i) = tilde_lopt(1,i)^2 * G(i,i) - 2*tilde_lopt(i)*z(i);
    end
    [min_cost, min_index]= min(cost);
    
    if min_cost < 0
        d(:,k+1) = d(:,k);
        d(min_index, k+1) = d_tilde(min_index,1);
        z = z - tilde_lopt(1,min_index)* G(:,min_index);
        k = k+1;
    else
        terminate = 1;
    end
end

%% --------------------------------------------------------
%%% Symbol to bit converion -------------------------------
%% --------------------------------------------------------
ML_symbol = d(:,k);
idxhat = common_find_idx(ML_symbol, par.RealSymbols);
% rearrange the detected symbols
idxhat1 = idxhat;
bithat = par.realbits(idxhat1,:);
bithat = reshape(bithat', par.Q, par.MT)';
end

