%% Hw #4 - Aeroelasticity - ME597/AAE556
% Victoria Nagorski - 12/2/22

%% Problem No. 1 Part a
clear; clc; close all;
syms h(t) theta(t) beta(t) t x_theta x_beta rho m x y I_theta I_beta c a...
    T1 T2 T3 T4 T5 T6 T7 T8 T9 T10 T11 T12 T13 T14 T15 T16 T17 T18 T19...
    C(k) b V Kh Kt Kb r_theta r_beta omega_h omega_t omega_b mu

variables = [omega_h omega_t  omega_b  a    c   b x_theta x_beta   r_theta   r_beta mu];
values =       [50    100      300   -0.4  0.6  1  0.2    0.0125  sqrt(0.25) sqrt(0.00625) 40]; 

B = 1;

% Define z
z = h + theta*x_theta + beta*x_beta;

% Give force and moment definitions
Q = V*theta + diff(h,t) + diff(theta,t)*b*((1/2)-a)+...
    (V/pi)*T10*beta + (b/(2*pi))*T11*diff(beta,t);
F = -pi*rho*b^2 *(diff(h,t,2)+V*diff(theta,t)-b*a*diff(theta,t,2)-...
    (V/pi)*T4*diff(beta,t)-(b/pi)*T1*diff(beta,t,2)) - 2*pi*rho*V*b*Q*C(k);
M = pi*rho*b^2*(b*a*diff(h,t,2)-V*b*((1/2)-a)*diff(theta,t)-...
    b^2*((1/8)+a^2)*diff(theta,t,2)-(V^2/pi)*(T4+T10)*beta+...
    ((V*b)/pi)*(-T1+T8+(c-a)*T4-(1/2)*T11)*diff(beta,t)+...
    b^2/pi*(T7+(c-a)*T1)*diff(beta,t,2)) + 2*pi*rho*V*b^2*(a+(1/2))*Q*C(k);
M_beta = pi*rho*b^2*((b/pi)*T1*diff(h,t,2)+((V*b)/pi)*(2*T9 +T1-...
    (a-(1/2))*T4)*diff(theta,t)-((2*b^2)/pi)*T13*diff(theta,t,2)-...
    (V/pi)^2*(T5-T4*T10)*beta+ ((V*b)/(2*pi*pi))*T4*T11*diff(beta,t)+...
    (b/pi)^2*T3*diff(beta,t,2)) - rho*V*b^2*T12*Q*C(k);

AeroTemp = expand([F/(m*b) M/(m*b^2) M_beta/(m*b^2)]);

% Solve for T (Kinetic Energy)
z_dot = diff(z,t);
temp = (1/2) * rho * int(int(expand(z_dot^2),x),y);
temp2 = expand(subs(expand(temp),rho*x*y,m));
temp3 = expand(subs(temp2,m*x_theta^2,m*r_theta^2));
T = subs(temp3,x_beta^2*m,m*r_beta^2);

% Solve for U (Potential Energy)
U = (1/2)*(Kh*h^2 + Kt*theta^2 + Kb*beta^2);

% Solve Lagranges Equations
mode = [h theta beta];
mode = mode(t);
eq = sym(zeros(length(mode),1));
for i = 1:length(mode)
    M = mode(i);
    dot = diff(M,t);
    eq(i,1) = simplify(diff(diff(T,dot),t) - diff(T,M) + diff(U,M));
end

damping = false;

% Create the matrices
sys.M = jacobian(eq,diff(mode,t,2)); % Mass
sys.K = jacobian(eq,mode); % 
if damping == false
    sys.C = zeros(length(mode),length(mode));
else
    sys.C = jacobian(eq,diff(mode,t));
end

% Plug in values for M and K (constants)
M = subs(sys.M,x_beta*x_theta,(r_beta^2 + x_beta*b*(c-a)));
M = double(subs(M/m,variables,values));
K = subs(sys.K,[Kh Kt Kb],[omega_h^2 omega_t^2*r_theta^2 omega_b^2*r_beta^2]);
K = double(subs(K,variables,values));

% Create the Aero matrices
Mu = m/(pi*rho*b^2);
Kk = b/V;
TempM = simplify(Mu*jacobian(AeroTemp,diff(mode,t,2))*diag([b 1 1]));
TempC = simplify(Kk*Mu*jacobian(AeroTemp,diff(mode,t))*diag([b 1 1]));
TempK = simplify(Kk^2*Mu*jacobian(AeroTemp,mode)*diag([b 1 1]));

% Definitions for the T values
c = 0.6;
a = -0.4;
t1 = -((2+c^2)/3)*sqrt(1-c^2)+c*acos(c);
t3 = -((1-c^2)/8)*(5*c^2+4)+(1/4)*c*(7+2*c^2)*sqrt(1-c^2)*acos(c)-((1/8)+c^2)*acos(c)^2;
t4 = c*sqrt(1-c^2)-acos(c);
t5 = -(1-c^2)-acos(c)^2+2*c*sqrt(1-c^2)*acos(c);
t7 = c*((7+2*c^2)/8)*sqrt(1-c^2)-((1/8)+c^2)*acos(c);
t8 = -(1/3)*(1+2*c^2)*sqrt(1-c^2)+c*acos(c);
t9 = (1/2)*((sqrt(1-c^2)*(1-c^2))/3 + a*t4);
t10 = sqrt(1-c^2) + acos(c);
t11 = (2-c)*sqrt(1-c^2)+(1-2*c)*acos(c);
t12 = (2+c)*sqrt(1-c^2)-(1+2*c)*acos(c);
t13 = -(1/2)*(t7+(c-a)*t1);
t15 = t4+t10;
t16 = t1-t8-(c-a)*t4+(1/2)*t11;
t17 = -2*t9-t1+(a-(1/2))*t4;
t18 = t5-t4*t10;
t19 = -(1/2)*t4*t11;
TNums = [T1 T3 T4 T5 T7 T8 T9 T10 T11 T12 T13 T15 T16 T17 T18 T19];
tNums = [t1 t3 t4 t5 t7 t8 t9 t10 t11 t12 t13 t15 t16 t17 t18 t19];

% Plug in the T values
TTempM = simplify(subs(TempM,TNums,tNums));
TTempC = simplify(subs(TempC,TNums,tNums));
TTempK = simplify(subs(TempK,TNums,tNums));

% Condense to one matrix before looping
AeroT = (-TTempM + (j/k)*TTempC + (1/(k^2))*TTempK)/mu;

% Begin Iteration
num = 700;
rst1 = zeros(3,num);
rst2 = zeros(3,num);
vel = zeros(3,num);
for int = num:-1:1
    rk = int*0.001;      % Iterate on the k value
    
    % Plug in values for aero
    if rk == 0
        CC = 1;
    else
        CC = besselk(1,j*rk)/(besselk(0,j*rk)+besselk(1,j*rk));
    end
    z = double(subs(AeroT,[variables k C(k)],[values rk CC]));
    f = real(z);
    g = imag(z);
    Aero = f + j*g;
    
    % Solve eigenvalue problem
    eigen = eig(inv(K) * (M + Aero));
    Reig = abs(real(eigen));
    Ieig = imag(eigen);
    rst1(:,int) = sqrt(1./Reig);
    rst2(:,int) = Ieig./Reig;
    vel(:,int) = sqrt(1./Reig)*B/rk;
end

% Graph results
figure
subplot(2,1,1)
hold on
plot(vel,rst1,'*r')
grid on
xlabel('V [m/s]')
ylabel('\omega')
axis([0 1000 0 1000])
hold off
subplot(2,1,2)
hold on
plot(vel,rst2,'*r')
grid on
xlabel('V [m/s]')
ylabel('g')             % Eq 4.8.13b in textbook
axis([0 500 -0.5 0.5])
line([0 500],[0 0],'Color','red','LineStyle','--')

%% Problem No. 1 Part b
clear; clc; close all;

% ======= Calc New Values ========
Vf = 301.5;
Vf10 = Vf*1.10;
Vf15 = Vf*1.15;
Vf20 = Vf*1.20;
delta_M = 8.75; % kg

omega_H = 50;
omega_T = 100;
omega_B = 300;
A = -0.4;
C= 0.6;
B = 1;
x_thetaB = 0.2;
x_betaB = 0.0125;
r_thetaB = sqrt(0.25);
r_betaB = sqrt(0.00625);
MU = 40;
RHO = 1.225;    % Density at Sea-Level [Kg/m3]

% Solve for the new c.g. w/point mass
m0 = MU*pi*RHO*B^2; % Solve for the initial mass
x_12 = B;       % Halfway mark
x_sc = x_12 - A*B;
x_cg0 = x_thetaB*B - A*B + x_12;
x2= 0.05*2*B;   % 5% of the cord (2*b)
x_cg_new = (m0*x_cg0 + x2*delta_M)/(m0 + delta_M);

% Calc the new x_theta bar
x_thetaB_new = (x_cg_new - x_sc)/B;

% Calc the new MU
MU_new = (m0 + delta_M)/(pi*RHO*B^2);

% Calc inertias and new r_theta bar
I_theta0 = m0*(B*r_thetaB)^2;
I0 = I_theta0 - m0*(x_thetaB*B)^2;
I_theta_new = I0 + (m0 + delta_M)*B^2*x_thetaB_new^2;
r_thetaB_new = sqrt(I_theta_new/(B^2*(m0 + delta_M)));

% Calc Kh and Ktheta (these do not change)
Kh = omega_H^2*m0;
Ktheta = omega_T^2*I_theta0;

% Calculate the new omegas
omega_H_new = sqrt(Kh/(m0 + delta_M));
omega_T_new = sqrt(Ktheta/I_theta_new);

% ========== Start V-g code Prep ============== 

syms h(t) theta(t) beta(t) t x_theta x_beta rho m x y I_theta I_beta c a...
    T1 T2 T3 T4 T5 T6 T7 T8 T9 T10 T11 T12 T13 T14 T15 T16 T17 T18 T19...
    C(k) b V Kh Kt Kb r_theta r_beta omega_h omega_t omega_b mu Mw ma

variables = [omega_h omega_t  omega_b  a    c   b x_theta x_beta   r_theta   r_beta mu Mw ma];
values = [omega_H_new  omega_T_new  300 -0.4  0.6  1  x_thetaB_new    0.0125  r_thetaB_new sqrt(0.00625) MU_new (m0+delta_M) m0]; 

% Define z
z = h + theta*x_theta + beta*x_beta;

% Give force and moment definitions
Q = V*theta + diff(h,t) + diff(theta,t)*b*((1/2)-a)+...
    (V/pi)*T10*beta + (b/(2*pi))*T11*diff(beta,t);
F = -pi*rho*b^2 *(diff(h,t,2)+V*diff(theta,t)-b*a*diff(theta,t,2)-...
    (V/pi)*T4*diff(beta,t)-(b/pi)*T1*diff(beta,t,2)) - 2*pi*rho*V*b*Q*C(k);
M = pi*rho*b^2*(b*a*diff(h,t,2)-V*b*((1/2)-a)*diff(theta,t)-...
    b^2*((1/8)+a^2)*diff(theta,t,2)-(V^2/pi)*(T4+T10)*beta+...
    ((V*b)/pi)*(-T1+T8+(c-a)*T4-(1/2)*T11)*diff(beta,t)+...
    b^2/pi*(T7+(c-a)*T1)*diff(beta,t,2)) + 2*pi*rho*V*b^2*(a+(1/2))*Q*C(k);
M_beta = pi*rho*b^2*((b/pi)*T1*diff(h,t,2)+((V*b)/pi)*(2*T9 +T1-...
    (a-(1/2))*T4)*diff(theta,t)-((2*b^2)/pi)*T13*diff(theta,t,2)-...
    (V/pi)^2*(T5-T4*T10)*beta+ ((V*b)/(2*pi*pi))*T4*T11*diff(beta,t)+...
    (b/pi)^2*T3*diff(beta,t,2)) - rho*V*b^2*T12*Q*C(k);

AeroTemp = expand([F/(Mw*b) M/(Mw*b^2) M_beta/(Mw*b^2)]);

% Solve for T (Kinetic Energy)
z_dot = diff(z,t);
temp = expand((1/2) * rho * z_dot^2);
temp2 = expand(subs(expand(temp),rho*x_theta^2,Mw*r_theta^2));
temp3 = expand(subs(temp2,rho*x_beta^2,ma*r_beta^2));
temp4 = expand(subs(temp3,rho*diff(h(t), t)^2,Mw*diff(h(t), t)^2));
temp5 = expand(subs(temp4,rho*x_beta*diff(beta(t), t)*diff(h(t), t),ma*x_beta*diff(beta(t), t)*diff(h(t), t)));
temp6 = expand(subs(temp5,rho*x_theta*diff(h(t), t)*diff(theta(t), t),Mw*x_theta*diff(h(t), t)*diff(theta(t),t)));
T = subs(temp6,rho*x_beta*x_theta*diff(beta(t), t)*diff(theta(t), t),(ma*r_beta^2 + ma*x_beta*(c-a))*diff(beta(t), t)*diff(theta(t), t));

% Solve for U (Potential Energy)
U = (1/2)*(Kh*h^2 + Kt*theta^2 + Kb*beta^2);

% Solve Lagranges Equations
mode = [h theta beta];
mode = mode(t);
eq = sym(zeros(length(mode),1));
for i = 1:length(mode)
    M = mode(i);
    dot = diff(M,t);
    eq(i,1) = simplify(diff(diff(T,dot),t) - diff(T,M) + diff(U,M));
end

damping = false;

% Create the matrices
sys.M = jacobian(eq,diff(mode,t,2)); % Mass
sys.K = jacobian(eq,mode); % 
if damping == false
    sys.C = zeros(length(mode),length(mode));
else
    sys.C = jacobian(eq,diff(mode,t));
end

% Plug in values for M and K (constants)
% M = subs(sys.M,x_beta*x_theta,(r_beta^2 + x_beta*b*(c-a)));
M = double(subs(expand(sys.M/Mw),variables,values));
K = subs(sys.K,[Kh Kt Kb],[omega_h^2 omega_t^2*r_theta^2 (omega_b^2*r_beta^2)/(1 + delta_M/m0)]);
K = double(subs(K,variables,values));

% Create the Aero matrices
Mu = Mw/(pi*rho*b^2);
Kk = b/V;
TempM = simplify(Mu*jacobian(AeroTemp,diff(mode,t,2))*diag([b 1 1]));
TempC = simplify(Kk*Mu*jacobian(AeroTemp,diff(mode,t))*diag([b 1 1]));
TempK = simplify(Kk^2*Mu*jacobian(AeroTemp,mode)*diag([b 1 1]));

% Definitions for the T values
c = 0.6;
a = -0.4;
t1 = -((2+c^2)/3)*sqrt(1-c^2)+c*acos(c);
t3 = -((1-c^2)/8)*(5*c^2+4)+(1/4)*c*(7+2*c^2)*sqrt(1-c^2)*acos(c)-((1/8)+c^2)*acos(c)^2;
t4 = c*sqrt(1-c^2)-acos(c);
t5 = -(1-c^2)-acos(c)^2+2*c*sqrt(1-c^2)*acos(c);
t7 = c*((7+2*c^2)/8)*sqrt(1-c^2)-((1/8)+c^2)*acos(c);
t8 = -(1/3)*(1+2*c^2)*sqrt(1-c^2)+c*acos(c);
t9 = (1/2)*((sqrt(1-c^2)*(1-c^2))/3 + a*t4);
t10 = sqrt(1-c^2) + acos(c);
t11 = (2-c)*sqrt(1-c^2)+(1-2*c)*acos(c);
t12 = (2+c)*sqrt(1-c^2)-(1+2*c)*acos(c);
t13 = -(1/2)*(t7+(c-a)*t1);
t15 = t4+t10;
t16 = t1-t8-(c-a)*t4+(1/2)*t11;
t17 = -2*t9-t1+(a-(1/2))*t4;
t18 = t5-t4*t10;
t19 = -(1/2)*t4*t11;
TNums = [T1 T3 T4 T5 T7 T8 T9 T10 T11 T12 T13 T15 T16 T17 T18 T19];
tNums = [t1 t3 t4 t5 t7 t8 t9 t10 t11 t12 t13 t15 t16 t17 t18 t19];

% Plug in the T values
TTempM = simplify(subs(TempM,TNums,tNums));
TTempC = simplify(subs(TempC,TNums,tNums));
TTempK = simplify(subs(TempK,TNums,tNums));

% Condense to one matrix before looping
AeroT = (-TTempM + (j/k)*TTempC + (1/(k^2))*TTempK)/mu;

% Begin Iteration
num = 700;
rst1 = zeros(3,num);
rst2 = zeros(3,num);
vel = zeros(3,num);
for int = num:-1:1
    rk = int*0.001;      % Iterate on the k value
    
    % Plug in values for aero
    if rk == 0
        CC = 1;
    else
        CC = besselk(1,j*rk)/(besselk(0,j*rk)+besselk(1,j*rk));
    end
    z = double(subs(AeroT,[variables k C(k)],[values rk CC]));
    f = real(z);
    g = imag(z);
    Aero = f + j*g;
    
    % Solve eigenvalue problem
    eigen = eig(inv(K) * (M + Aero));
    Reig = abs(real(eigen));
    Ieig = imag(eigen);
    rst1(:,int) = sqrt(1./Reig);
    rst2(:,int) = Ieig./Reig;
    vel(:,int) = sqrt(1./Reig)*B/rk;
end

% Graph results
figure
subplot(2,1,1)
hold on
plot(vel,rst1,'*r')
grid on
xlabel('V [m/s]')
ylabel('\omega')
axis([0 1000 0 1000])
hold off
subplot(2,1,2)
hold on
plot(vel,rst2,'*r')
grid on
xlabel('V [m/s]')
ylabel('g')             % Eq 4.8.13b in textbook
axis([300 370 -0.5 0.5])
line([0 500],[0 0],'Color','red','LineStyle','--')

figure
hold on
plot([5.25 7.15 8.76],[331.65 346.73 361.80])
xlabel('Delta M [kg]')
ylabel('V_F [m/s]')
hold off

%% Problem No. 2
clear;clc;close all;

% =================== Set-up the Problem ========================
syms h(t) theta(t) beta(t) t x_theta x_beta rho m x y I_theta I_beta c a...
    T1 T2 T3 T4 T5 T6 T7 T8 T9 T10 T11 T12 T13 T14 T15 T16 T17 T18 T19...
    C(k) b V Kh Kt Kb r_theta r_beta omega k omega_h omega_t mu Lh La ...
    Ma d_bar Cla

% Define z
z = h + theta*x;

% Solve for T (Kinetic Energy)
z_dot = diff(z,t);
temp = expand((1/2) * rho * z_dot^2);
temp2 = expand(subs(temp,rho*x^2,m*r_theta^2));
temp3 = expand(subs(temp2,rho*x,m*x_theta));
T = expand(subs(temp3,rho,m));

% Define U
U = (1/2)*Kh*h^2 + (1/2)*Kt*theta^2;

% Define the modes (h is really h/b)
mode = [h theta];

% Solve Lagrange's Equation
mode = mode(t);
eq = sym(zeros(length(mode),1));
for i = 1:length(mode)
    M = mode(i);
    dot = diff(M,t);
    eq(i,1) = simplify(diff(diff(T,dot),t) - diff(T,M) + diff(U,M));
end

% Change out defined terms in equations
eq2 = expand(subs(expand(eq/m),Kh/m,omega_h^2));
eq3 = expand(subs(eq2,Kt/m,omega_t^2*r_theta^2));

damping = false;

% Create the matrices
sys.M = jacobian(eq3,diff(mode,t,2)); % Mass
sys.K = jacobian(eq3,mode); % 
if damping == false
    sys.C = zeros(length(mode),length(mode));
else
    sys.C = jacobian(eq3,diff(mode,t));
end

% Solve for the aero matrix (h is h/b and theta)
Mh = 1/2;
P_mb = (1/mu)* (Lh*h + (La - (1/2 + a)*Lh)*theta);
M_mbb = (1/mu)* ((Mh-(1/2 + a)*Lh)*h + (Ma -(1/2 + a)*(La+Mh)+(1/2+a)^2*Lh)*theta);
aero_eq = [P_mb;M_mbb];
sys.Aero = jacobian(aero_eq,mode);

% ======================= The p-k Code ==========================
% Parameters
variables = [omega_h omega_t  a   d_bar b x_theta r_theta Cla mu];
values =       [10      25  -0.2   0.4  3  0.1  sqrt(0.5) 2*pi 20]; 
MU = 20;
B = 3;
OMEGAT = 25;

% Plug Into Mass and K matrices (Constant)
M = double(subs(sys.M,variables,values));
KK = double(subs(sys.K,variables,values));

% Replace Lh, La, and Ma in Aero pg. 415 and 417
lh = 1 - j*2*C(k)/k;
la = (1/2) - j*(1+2*C(k))/k - (2*C(k)/k^2);
ma = (3/8)- j/k;
AeroT = subs(mu*sys.Aero,[Lh La Ma],[lh la ma]);

% Begin Iteration Code
Nv = 1000;
Vel = linspace(0.001,1000,Nv);
sigma1 = zeros(1,Nv);
w1 = zeros(1,Nv);
sigma2 = zeros(1,Nv);
w2 = zeros(1,Nv);
for iter = 1:Nv
    V = Vel(iter);
    K = 0.5;            % Initial guess for k
    e = 1;              % Initialize error for k
    while e > 0.0001
        Ck = besselh(1,2,K)/(besselh(1,2,K)+1i*besselh(0,2,K));
        
        % Plug in values to Aero matrix
        Aero = double(subs(AeroT,[C(k) mu a k],[Ck 20 -0.2 K]));
        % Break Aero into Aij and Bij matricies
        Aij = V^2*K^2/(MU*B^2) * real(Aero);
        Bij = V/B *(K/MU)*imag(Aero);
        % Combine Matricies
        Kij_bar =  KK- Aij;

        % Make the Plant
        A = [zeros(length(M)) eye(length(M));
             -inv(M)*Kij_bar   inv(M)*Bij];

        % Solve eigenvalues
        p = eig(A);

        [wu index] = max(abs(imag(p)));
        if wu == 0
            e = 1e-6;
        else
            knew = wu*B/V;
            e = abs(knew-K);
            K = knew;
        end
    end
    sigma1(:,iter) = real(p(index))/imag(p(index));
    w1(:,iter) = imag(p(index));


    e = 1;              % Initialize error for force
    K = 0.5;            % Initilaize k
    while e > 0.0001
        Ck = besselh(1,2,K)/(besselh(1,2,K)+1i*besselh(0,2,K));
        
        % Plug in values to Aero matrix
        Aero = double(subs(AeroT,[C(k) mu a k],[Ck 20 -0.2 K]));
        % Break Aero into Aij and Bij matricies
        Aij = V^2*K^2/(MU*B^2) * real(Aero);
        Bij = V/B *(K/MU)*imag(Aero);
        % Combine Matricies
        Kij_bar = KK - Aij;

        % Make the Plant
        A = [zeros(length(M)) eye(length(M));
             -inv(M)*Kij_bar   inv(M)*Bij];

        % Solve eigenvalues
        p = eig(A);

        [wu1 index] = min(abs(imag(p)));
        if wu1 == 0
            e = 1e-6;
        else
            knew = wu1*B/V;
            e = abs(knew-K);
            K = knew;
        end
    end
    sigma2(:,iter) = real(p(index))/imag(p(index));
    w2(:,iter) = imag(p(index));
end

figure
subplot(2,1,1)
hold on
title('p-k Method')
plot(Vel,w1./OMEGAT,'r');
plot(Vel,w2./OMEGAT,'b')
grid on
xlabel('V [m/s]')
ylabel('\omega')
axis([0 250 0 1.25])
hold off
subplot(2,1,2)
hold on
plot(Vel,sigma1,'r');
plot(Vel,sigma2,'b')
grid on
xlabel('V [m/s]')
ylabel('g')             % Eq 4.8.13b in textbook
axis([0 250 -0.5 0.15])
line([0 250],[0 0],'Color','red','LineStyle','--')
hold off

% ======================= The v-g Code ==========================
% Plug in values for M and K (constants)
K = KK;

% Begin Iteration
num = 400;
rst1 = zeros(2,num);
rst2 = zeros(2,num);
vel = zeros(2,num);
for int = num:-1:1
    rk = int*0.01;      % Iterate on the k value
    
    % Plug in values for aero
    if rk == 0
        CC = 1;
    else
        CC = besselk(1,j*rk)/(besselk(0,j*rk)+besselk(1,j*rk));
    end
    z = double(subs(AeroT/MU,[variables k C(k)],[values rk CC]));
    f = real(z);
    g = imag(z);
    Aero = f + j*g;
    
    % Solve eigenvalue problem
    eigen = eig(inv(K) * (M + Aero));
    Reig = abs(real(eigen));
    Ieig = imag(eigen);
    rst1(:,int) = sqrt(1./Reig);
    rst2(:,int) = Ieig./Reig;
    vel(:,int) = sqrt(1./Reig)*B/rk;
end

% Graph results
figure
subplot(2,1,1)
hold on
title('V-g Method')
plot(vel,rst1,'*r')
grid on
xlabel('V [m/s]')
ylabel('\omega')
axis([0 250 0 50])
hold off
subplot(2,1,2)
hold on
plot(vel,rst2,'*r')
grid on
xlabel('V [m/s]')
ylabel('g')             % Eq 4.8.13b in textbook
axis([0 250 -0.5 0.15])
line([0 250],[0 0],'Color','red','LineStyle','--')

% ======================= Quasi-Steady ==========================
syms q V e
% Create the matrices
% e = 0.4;
rho = 1.225;    % Density at Sea-Level [Kg/m3]
LIFTq = q*Cla*b*(diff(h,t)/V + b/V*(1/2-a)*diff(theta,t) + theta);
aero_q = expand([LIFTq;-e*LIFTq]);

K = sys.K/omega_t^2;

damping = false;    % Set C=0 so that motion is harmonic (pg 348)
AA = jacobian(aero_q,mode);
AA = AA*(Cla*V^2)/(Cla*q*b*b^2*omega_t^2*pi*mu);
if damping == false
    sys.C = zeros(length(mode),length(mode));
else
    sys.C = jacobian(aero_q,diff(mode,t));
end
A_quasi = [zeros(length(mode),length(mode))    eye(length(mode));
                -inv(sys.M)*(K+AA)         -inv(sys.M)*sys.C];

A_quasi = simplify(subs(A_quasi,q,(1/2)*rho*V^2));
A_quasi = subs(A_quasi,[variables e],[values 0.3]);


% Plug in velocities
Velocities = 0:0.1:300;
other.frequencies = zeros(length(mode)*2,length(Velocities));
other.damping = zeros(length(mode)*2,length(Velocities));
for i = 1:length(Velocities)
    v = Velocities(i);
    Temp = double(subs(A_quasi,V,v));
    p = eig(Temp);
    other.frequencies(:,i) = imag(p);
    other.damping(:,i) = real(p);
end

% Plot everything
figure
hold on
for i = 1:(length(mode)*2)
     plot(Velocities,other.frequencies(i,:))
end
title('Frequency vs Flow Speed (Quasi-Steady)')
xlabel('Flow Speed')
ylabel('Frequency')
grid('on')
hold off

figure
hold on
for i = 1:(length(mode)*2)
     plot(Velocities,other.damping(i,:))
end
title('Damping vs Flow Speed (Quasi-Steady)')
xlabel('Flow Speed')
ylabel('Damping')
grid('on')
hold off