% number of flops vs number of users plot

eta = 4;
NT=  8:2:64;
 NR =4*NT;
% NR =128;
detector = {'RI','Cheby', 'GS', 'NS-RI', 'NS-GS', 'NS-Cheby',...
             'ST-Jacobi', 'CG'};


mult(1,:) = 28*NT.*NR + 16*NT;                 % RI
mult(2,:) = 36*NT.*NR + 48*NT;              % Cheby
mult(3,:) = 4*NT.*NR + 4*NT +12*NT.^2;         % GS
mult(4,:) = 44*NT.*NR + 24*NT;            % NS-RI
mult(5,:) = 28*NT.*NR + 8*NT.^2 + 16*NT;      % NS-GS
mult(6,:) = 52*NT.*NR + 48*NT;            % NS-Cheby
mult(7,:) = 20*NT.*NR + 8*NT.^2 + 32*NT;  % St jacobi
mult(8,:) = 28*NT.*NR + 84*NT;                  % CG







marker_style = {'bo-','rs--','mv-.','kp:','g*-','c>--','yx:', 'kx-',...
                'kh:', 'b+--', 'rs-', 'gd:', 'ys--', 'c+--', 'k*', 'kh:'};
figure(1)
for d=1:length(detector)
    if d==1
        plot(NT,mult(d,:),marker_style{d},'LineWidth',2)
        hold on
    else
        
        plot(NT,mult(d,:),marker_style{d},'LineWidth',2)
    end
end

hold off
grid on
xlabel('NT','FontSize',12)
ylabel('Number of numtiplication','FontSize',12)
% axis([min(par.SNRdB_list) max(par.SNRdB_list) 1e-6 1])
legend(detector,'FontSize',12)
set(gca,'FontSize',12)



%% iteration 4

mult1(1,:) = 36*NT.*NR + 20*NT;                 % RI
mult1(2,:) = 44*NT.*NR + 60*NT;              % Cheby
mult1(3,:) = 4*NT.*NR + 4*NT +16*NT.^2;         % GS
mult1(4,:) = mult(4,:);            % NS-RI
mult1(5,:) = mult(5,:);      % NS-GS
mult1(6,:) = mult(6,:);            % NS-Cheby
mult1(7,:) = 20*NT.*NR + 12*NT.^2 + 36*NT;  % St jacobi
mult1(8,:) = 36*NT.*NR + 112*NT;                  % CG
        

marker_style = {'bo-','rs--','mv-.','kp:','g*-','c>--','yx:', 'kx-',...
                'kh:', 'b+--', 'rs-', 'gd:', 'ys--', 'c+--', 'k*', 'kh:'};
figure(2)
for d=1:length(detector)
    if d==1
        plot(NT,mult1(d,:),marker_style{d},'LineWidth',2)
        hold on
    else
        
        plot(NT,mult1(d,:),marker_style{d},'LineWidth',2)
    end
end

hold off
grid on
xlabel('NT','FontSize',12)
ylabel('Number of numtiplication','FontSize',12)
% axis([min(par.SNRdB_list) max(par.SNRdB_list) 1e-6 1])
legend(detector,'FontSize',12)
set(gca,'FontSize',12)






NT=  8:2:64;

NT = NT';
 NR =4*NT;
 
RI = 28*NT.*NR + 16*NT;                 % RI
Cheby = 36*NT.*NR + 48*NT;              % Cheby
GS = 4*NT.*NR + 4*NT +12*NT.^2;         % GS
NS_RI = 44*NT.*NR + 24*NT;            % NS-RI
NS_GS = 28*NT.*NR + 8*NT.^2 + 16*NT;      % NS-GS
NS_Cheby = 52*NT.*NR + 48*NT;            % NS-Cheby
ST_Jacobi = 20*NT.*NR + 8*NT.^2 + 32*NT;  % St jacobi
CG = 28*NT.*NR + 84*NT;                  % CG

 
FileName = 'Compiled_Data_Origin/complexity.mat';
save(FileName,'NT','RI', 'Cheby',...
              'GS', 'NS_RI',... 
              'NS_GS', 'NS_Cheby','ST_Jacobi','CG');
