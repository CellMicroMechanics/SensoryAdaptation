clc;clear;

%% Varying Chemical Potential
% dt = 0.0001; 
% T = 20; % Simulation Time
% sysparam.dt = dt;
% N = T/dt+1;
% sysparam.N = N;
% 
% omega_m = 5;
% omega_a = 50;
% sigma_m = 0.1;
% sigma_a = 0.1;
% delta_a = 0.01*omega_a;
% delta_m = 0.01*omega_m;
% sysparam.omega_m = omega_m;
% sysparam.omega_a = omega_a;
% sysparam.delta_a = delta_a;
% sysparam.delta_m = delta_m;
% 
% 
% C = delta_m*omega_a/omega_m/delta_a;
% a0 = 0.5;
% m = 4;
% H = 1;
% 
% sysparam.H = H;
% sysparam.mmax = 4;
% sysparam.a0 = a0;
% temp = delta_m/omega_m;
% sysparam.C = C;
% 
% Beta = [0:0.02:1];
% s = 20;
% bb = length(Beta);
% G = @(m) 1./(1 + (s./exp(2*m)).^H);
% Fa = @(a,G) -omega_a.*(a - G);
% dG = @(s,m) (2*exp(-2*m)*H.*s.*(exp(-2*m).*s).^(-1+H))./(1+(s.*exp(-2*m)).^H).^2;
% 
% n = 1e4;
% sysparam.s = s;
% m_str = -1/2*log(1/s*(1/a0-1));
% 
% for k = 1:bb
%     tic
%     beta = Beta(k)
%     sysparam.beta = beta;
%     Fm = @(a0,a,dG) -omega_m.*(a-a0).*(beta - (1-beta).*C.*dG);
% 
%     l = N+1;
%     
%     sum_Qa = zeros([1,l-2]);
%     sum_Qm = zeros([1,l-2]);
%     Fa_str = zeros([1,l-2]);
%     Fm_str = zeros([1,l-2]); % Stratonovich Forces
%     A_mean = 0;
%     M_mean = 0;
%     parfor j = 1:n
%        
%         [A,M] = gen_traj(sysparam);
%         
%         A_mean = A_mean + mean(A(end-1e5:end))./n;
%          M_mean = M_mean + mean(M(end-1e5:end))./n;
%         % Energetic Computations 
%         i = 2:l-1;
%         a_str = (A(i)+A(i-1))/2;
%         m_str = (M(i)+M(i-1))/2;
% 
%         Fa_str = Fa(a_str,G(m_str));
%         Fm_str = Fm(a0, a_str, dG(s,m_str));
% 
%         dadt = (A(i)-A(i-1))/dt ;
%         dmdt = (M(i)-M(i-1))/dt;
% 
%         sum_Qa = sum_Qa + Fa_str.*dadt/delta_a;
%         sum_Qm = sum_Qm + Fm_str.*dmdt/delta_m;
%         
%         if (mod(j,100) == 0)
%             j
%         end
%         
%     end
%     m_star(k) = M_mean;
%     Err(k) = abs(A_mean-a0)
%     
%     AvgQa0 = sum_Qa./n;
%     AvgQm0 = sum_Qm./n;
%     Diss(k) = omega_m*beta^2*temp*(1+dG(s,m_star(k)))^2
%     diss(k) = temp*mean(AvgQa0(end-1e5:end)+ AvgQm0(end-1e5:end))
%     
%     toc
% end
% 
% beta_c = (C*dG(s,m_str))/(C*dG(s,m_str) + 1)
% 
% figure
% plot(Diss,Err,'-o')
% hold on
% plot(diss,Err,'--')
% grid on
% legend('$\dot{W}^a_{diss}$','$\dot{W}_{diss}$','Interpreter','latex')


%% Varying input

% clc;clear;
% dt = 0.0001; 
% T = 20; % Simulation Time
% sysparam.dt = dt;
% N = T/dt+1;
% sysparam.N = N;
% beta = 1;
% 
% omega_m = 5;
% omega_a = 50;
% sigma_m = 0.1;
% sigma_a = 0.1;
% delta_a = 0.01*omega_a;
% delta_m = 0.01*omega_m;
% sysparam.omega_m = omega_m;
% sysparam.omega_a = omega_a;
% sysparam.delta_a = delta_a;
% sysparam.delta_m = delta_m;
% sysparam.beta = beta;
% 
% 
% C = delta_m*omega_a/omega_m/delta_a;
% a0 = 0.5;
% m = 4;
% H = 1;
% 
% sysparam.H = H;
% sysparam.mmax = 4;
% sysparam.a0 = a0;
% temp = delta_m/omega_m;
% sysparam.C = C;
% 
% S = [250:250:2e4]
% ss = length(S);
% 
% n = 1e4;
% 
% for k = 1:ss
%     tic
%     s = S(k)
%     sysparam.s = s;
%     G = @(m) 1./(1 + (s./exp(2*m)).^H);
%     Fa = @(a,G) -omega_a.*(a - G);
%     dG = @(s,m) (2*exp(-2*m)*H.*s.*(exp(-2*m).*s).^(-1+H))./(1+(s.*exp(-2*m)).^H).^2;
%     Fm = @(a0,a,dG) -omega_m.*(a-a0).*(beta - (1-beta).*C.*dG);
% 
%     l = N+1;
%     
%     sum_Qa = zeros([1,l-2]);
%     sum_Qm = zeros([1,l-2]);
%     Fa_str = zeros([1,l-2]);
%     Fm_str = zeros([1,l-2]); % Stratonovich Forces
%     A_mean = 0;
%     M_mean = 0;
%     parfor j = 1:n
%        
%         [A,M] = gen_traj(sysparam);
%         
%         A_mean = A_mean + mean(A(end-1e5:end))./n;
%          M_mean = M_mean + mean(M(end-1e5:end))./n;
%         % Energetic Computations 
%         i = 2:l-1;
%         a_str = (A(i)+A(i-1))/2;
%         m_str = (M(i)+M(i-1))/2;
% 
%         Fa_str = Fa(a_str,G(m_str));
%         Fm_str = Fm(a0, a_str, dG(s,m_str));
% 
%         dadt = (A(i)-A(i-1))/dt ;
%         dmdt = (M(i)-M(i-1))/dt;
% 
%         sum_Qa = sum_Qa + Fa_str.*dadt/delta_a;
%         sum_Qm = sum_Qm + Fm_str.*dmdt/delta_m;
%         
%         if (mod(j,100) == 0)
%             j
%         end
%         
%     end
%     m_star(k) = M_mean;
%     Err(k) = abs(A_mean-a0)
%     
%     AvgQa0 = sum_Qa./n;
%     AvgQm0 = sum_Qm./n;
%     Diss(k) = omega_m*beta^2*temp*(1+dG(s,m_star(k)))^2
%     diss(k) = temp*mean(AvgQa0(end-1e5:end)+ AvgQm0(end-1e5:end))
%     m_str(k) = -1/2*log(1/s*(1/a0-1));
%     toc
% end
% 
% 
% figure
% plot(Diss,Err,'-o')
% hold on
% plot(diss,Err,'--')
% grid on
% legend('$\dot{W}^a_{diss}$','$\dot{W}_{diss}$','Interpreter','latex')
