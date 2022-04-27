function [P, norm2] = generate_pilot(t,N)

%% DFT based Pilot generation t x N

Wn = exp(1i*2*pi/N);

P = zeros(t,N);

for i  = 1:t
    for j=1:N
        
    P(i,j) = Wn^((i-1)*(j-1))  ;  
    end
end

norm2 = norm(P)^2;
P = sqrt(norm2/(t*N))*P;
% pilot_power_dft = norm(P(:,2))^2
%% random orthogonal Pilot generation t x N
% pilot = (1/sqrt(2))*(randn(N,t) + 1i*randn(N,t));
% [Q, R] = qr(pilot);
% Xp_tmp = Q(:,1:t).';
% pilot_power_qr = norm(Xp_tmp(:,2))^2
%Xp_tmp*Xp_tmp'

