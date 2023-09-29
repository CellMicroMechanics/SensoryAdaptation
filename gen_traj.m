function [A,M] = gen_traj(sysparam)

dt = sysparam.dt;
omega_a = sysparam.omega_a;
omega_m = sysparam.omega_m;
delta_a = sysparam.delta_a;
delta_m = sysparam.delta_m;
C = sysparam.C;
H = sysparam.H;
beta = sysparam.beta;
N = sysparam.N;
l = N+1;
s = sysparam.s;
a0 = sysparam.a0;
mmax = sysparam.mmax;
%% Functions
G = @(m) 1./(1 + (s./exp(2*m)).^H);
Fa = @(a,G) -omega_a.*(a - G);
dG = @(m) (2*exp(-2*m)*H.*s.*(exp(-2*m).*s).^(-1+H))./(1+(s.*exp(-2*m)).^H).^2;
Fm = @(a0,a,dG) -omega_m.*(a-a0).*(beta - (1-beta).*C.*dG);
%%

A_init = ones([1,l]).*0.5; % Initial Condition a
M_init = ones([1,l]).*0.5; % Initial Condition m
A = A_init; % Initialize A
M = M_init; % Initialize M
FA = zeros([1,l-1]);
FM = zeros([1,l-1]);
pdf= makedist('Normal',0,sqrt(dt));
for i = 2:l
    
   dWt_a = sqrt(2*delta_a)*random(pdf);
   dWt_m = sqrt(2*delta_m)*random(pdf);

   FA(i-1) = Fa(A(i-1), G(M(i-1)));
   FM(i-1) = Fm(a0,A(i-1), dG(M(i-1)));

   A(i) = A(i-1) + FA(i-1)*dt + dWt_a;
   M(i) = M(i-1) + FM(i-1)*dt + dWt_m;


    %Boundary Condition on m
    if (abs(M(i)) <= mmax)
        M(i) = abs(M(i));
    else
        M(i) = mod(M(i),2*mmax); % abs of central modulus
        if M(i) > mmax
           M(i) = abs(M(i) - 2*mmax);
        end % or simply use abs(modcent(X(2,i),2*m))
    end
    % Boundary Condition on a
    if (abs(A(i)) <= 1)
        A(i) = abs(A(i));
    else 
        A(i) = mod(A(i),2);
        if A(i) > 1
           A(i) = abs(A(i) - 2);
        end
    end

end

end