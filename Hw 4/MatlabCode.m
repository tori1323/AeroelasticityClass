%% Hw #4 - Aeroelasticity - ME597/AAE556
% Victoria Nagorski - 11/13/22


%% Problem No.1
clear;clc;close all

% Define the variables
b_c = 6;
e_c = 0.1;
c_e = 1/e_c;
b_e = b_c/e_c;
c = 2;  % Barroed from problem 2
e = e_c*c;
b = b_c*c;
E = 0.15;
K_phi_theta = 3;
V = 50;
rho = 1.225;
a0 = 2*pi;
Q = (1/2)*rho*V^2*c*b*a0;
K_theta = 1600;     % Found this value online
Q_ratio = 0.5;

% Solve for coeffiecient ratios with Glauert's theory
c_l_a = (1/pi) * (acos(1-(2*E)) + 2*sqrt(E*(1-E)));
c_m_a = (-1/pi) * (1-E) * sqrt((1-E)*E);

% Variables to varry in graphs
Lambda = [-30 0 30].*pi/180;
gamma = (-30:0.01:30).*pi/180;

% Define the equation
bottom = K_theta.*K_phi_theta.*(cos(gamma).^2 + sin(gamma).^2).^2.*(1-Q_ratio);
L_F_R = @(L) Q*e*[(K_phi_theta*cos(gamma).^2+sin(gamma).^2+...
    (K_phi_theta-1).*sin(gamma).*cos(gamma).*tan(L)).*(1 + ...
    c_m_a/c_l_a*c_e) - ((K_phi_theta-1).*sin(gamma).*cos(gamma)+...
    K_phi_theta.*sin(gamma).^2.*tan(L)+...
    cos(gamma).^2*tan(L)).*(b_e/2)]./bottom + 1;

% Graphing
figure
hold on
for i = 1:length(Lambda)
    lambda = Lambda(1,i);          % Pull out the q_hat value being used
    temp_ratio = L_F_R(lambda);        % Solve for the K_e values
    plot(gamma.*180/pi,temp_ratio,'LineWidth',3)         % Graph everything
end
title('$\frac{L_F}{L_R}$ vs $\gamma$ ($\frac{Q}{Q_D}$ = 0.5)','Interpreter','latex')
xlabel('$\gamma$','Interpreter','latex')
ylabel('$\frac{L_F}{L_R}$','Interpreter','latex')
legend('$\Lambda$ = $-30^o$','$\Lambda$ = $0.0^o$','$\Lambda$ = $30^o$','Interpreter','latex')
grid('on')
hold off

Q_ratio = 0.75;

% Solve for coeffiecient ratios with Glauert's theory
c_l_a = (1/pi) * (acos(1-(2*E)) + 2*sqrt(E*(1-E)));
c_m_a = (-1/pi) * (1-E) * sqrt((1-E)*E);

% Variables to varry in graphs
Lambda = [-30 0 30].*pi/180;
gamma = (-30:0.01:30).*pi/180;

% Define the equation
bottom = K_theta.*K_phi_theta.*(cos(gamma).^2 + sin(gamma).^2).^2.*(1-Q_ratio);
L_F_R = @(L) Q*e*[(K_phi_theta*cos(gamma).^2+sin(gamma).^2+...
    (K_phi_theta-1).*sin(gamma).*cos(gamma).*tan(L)).*(1 + ...
    c_m_a/c_l_a*c_e) - ((K_phi_theta-1).*sin(gamma).*cos(gamma)+...
    K_phi_theta.*sin(gamma).^2.*tan(L)+...
    cos(gamma).^2*tan(L)).*(b_e/2)]./bottom + 1;

% Graphing
figure
hold on
for i = 1:length(Lambda)
    lambda = Lambda(1,i);          % Pull out the q_hat value being used
    temp_ratio = L_F_R(lambda);        % Solve for the K_e values
    plot(gamma.*180/pi,temp_ratio,'LineWidth',3)         % Graph everything
end
title('$\frac{L_F}{L_R}$ vs $\gamma$ ($\frac{Q}{Q_D}$ = 0.75)','Interpreter','latex')
xlabel('$\gamma$','Interpreter','latex')
ylabel('$\frac{L_F}{L_R}$','Interpreter','latex')
legend('$\Lambda$ = $-30^o$','$\Lambda$ = $0.0^o$','$\Lambda$ = $30^o$','Interpreter','latex')
grid('on')
hold off

Q_ratio = 1.0;

% Solve for coeffiecient ratios with Glauert's theory
c_l_a = (1/pi) * (acos(1-(2*E)) + 2*sqrt(E*(1-E)));
c_m_a = (-1/pi) * (1-E) * sqrt((1-E)*E);

% Variables to varry in graphs
Lambda = [-30 0 30].*pi/180;
gamma = (-30:0.01:30).*pi/180;

% Define the equation
bottom = K_theta.*K_phi_theta.*(cos(gamma).^2 + sin(gamma).^2).^2.*(1-Q_ratio);
L_F_R = @(L) Q*e*[(K_phi_theta*cos(gamma).^2+sin(gamma).^2+...
    (K_phi_theta-1).*sin(gamma).*cos(gamma).*tan(L)).*(1 + ...
    c_m_a/c_l_a*c_e) - ((K_phi_theta-1).*sin(gamma).*cos(gamma)+...
    K_phi_theta.*sin(gamma).^2.*tan(L)+...
    cos(gamma).^2*tan(L)).*(b_e/2)]./bottom + 1;

% Graphing
figure
hold on
for i = 1:length(Lambda)
    lambda = Lambda(1,i);          % Pull out the q_hat value being used
    temp_ratio = L_F_R(lambda);        % Solve for the K_e values
    plot(gamma.*180/pi,temp_ratio,'LineWidth',3)         % Graph everything
end
title('$\frac{L_F}{L_R}$ vs $\gamma$ ($\frac{Q}{Q_D}$ = 1.0)','Interpreter','latex')
xlabel('$\gamma$','Interpreter','latex')
ylabel('$\frac{L_F}{L_R}$','Interpreter','latex')
legend('$\Lambda$ = $-30^o$','$\Lambda$ = $0.0^o$','$\Lambda$ = $30^o$','Interpreter','latex')
grid('on')
hold off

Q_ratio = 1.25;

% Solve for coeffiecient ratios with Glauert's theory
c_l_a = (1/pi) * (acos(1-(2*E)) + 2*sqrt(E*(1-E)));
c_m_a = (-1/pi) * (1-E) * sqrt((1-E)*E);

% Variables to varry in graphs
Lambda = [-30 0 30].*pi/180;
gamma = (-30:0.01:30).*pi/180;

% Define the equation
bottom = K_theta.*K_phi_theta.*(cos(gamma).^2 + sin(gamma).^2).^2.*(1-Q_ratio);
L_F_R = @(L) Q*e*[(K_phi_theta*cos(gamma).^2+sin(gamma).^2+...
    (K_phi_theta-1).*sin(gamma).*cos(gamma).*tan(L)).*(1 + ...
    c_m_a/c_l_a*c_e) - ((K_phi_theta-1).*sin(gamma).*cos(gamma)+...
    K_phi_theta.*sin(gamma).^2.*tan(L)+...
    cos(gamma).^2*tan(L)).*(b_e/2)]./bottom + 1;

% Graphing
figure
hold on
for i = 1:length(Lambda)
    lambda = Lambda(1,i);          % Pull out the q_hat value being used
    temp_ratio = L_F_R(lambda);        % Solve for the K_e values
    plot(gamma.*180/pi,temp_ratio,'LineWidth',3)         % Graph everything
end
title('$\frac{L_F}{L_R}$ vs $\gamma$ ($\frac{Q}{Q_D}$ = 1.25)','Interpreter','latex')
xlabel('$\gamma$','Interpreter','latex')
ylabel('$\frac{L_F}{L_R}$','Interpreter','latex')
legend('$\Lambda$ = $-30^o$','$\Lambda$ = $0.0^o$','$\Lambda$ = $30^o$','Interpreter','latex')
grid('on')
hold off

%% Problem No.2 Part a
clear;clc;close all

% Define the variables
s = 7.5;                % [m]
c = 2;                  % [m]
xf = 0.48*c;            % [m]
m = 200;                % [kg/m^2]
EI = 2e7;               % [Nm^2]
GJ = 2e6;               % [Nm^2]
aw = 2*pi;              % Lift curve
rho = 1.225;            % [kg/m^3]
ab = xf - c/2;
e = (c/4 + ab)/c;
M_thetad = -1.2;         % Non-dimensional pitch damping derivative

syms V lamda
q = (1/2)*rho*V^2;      % Dynamic pressure

% Set up matrices (Mx_ddot + Dx_dot + Kx = 0)
M = m.*[(s/5)*c                     (s/4)*(((c^2)/2) - xf*c);
        (s/4)*(((c^2)/2) - xf*c)    (s/3)*(((c^3)/3) - xf*(c^2) + (xf^2)*c)];

D = [(1/5)*q*c*aw*s/V                       0;
     (-1/4)*q*c^2*e*aw*s/V      (-s/12)*q*c^3*M_thetad/V];

K = [4*EI/s^3       q*c*aw*(s/4);
        0           GJ/s - q*c^2*e*aw*(s/3)];

% Reduce order of equations
A = [zeros(2,2)   eye(2);
     -inv(M)*K  -inv(M)*D];

% Plug in velocities
Velocities = 0:1:180;
Frequencies = zeros(4,length(Velocities));
Damping = zeros(4,length(Velocities));
for i = 1:length(Velocities)
    v = Velocities(i);
    Temp = double(subs(A,V,v));
    [wn,zeta] = damp(Temp);
    Frequencies(:,i) = wn;
    Damping(:,i) = zeta;
end

% Plot everything for w/Aerodynamic Damping
figure
hold on
plot(Velocities,Frequencies(1,:))
plot(Velocities,Frequencies(2,:))
plot(Velocities,Frequencies(3,:))
plot(Velocities,Frequencies(4,:))
title('Frequency vs Flow Speed (w/Aerodynamic Damping)')
xlabel('Flow Speed')
ylabel('Frequency')
grid('on')
hold off

figure
hold on
plot(Velocities,100.*Damping(1,:))
plot(Velocities,100.*Damping(2,:))
plot(Velocities,100.*Damping(3,:))
plot(Velocities,100.*Damping(4,:))
title('Damping vs Flow Speed (w/Aerodynamic Damping)')
xlabel('Flow Speed')
ylabel('Damping (%)')
grid('on')
hold off

% Reduce order of equations
A = [zeros(2,2) eye(2);
     -inv(M)*K  zeros(2,2)];

% Plug in velocities
Velocities = 0:1:180;
Frequencies = zeros(4,length(Velocities));
Damping = zeros(4,length(Velocities));
for i = 1:length(Velocities)
    v = Velocities(i);
    Temp = double(subs(A,V,v));
    [wn,zeta] = damp(Temp);
    Frequencies(:,i) = wn;
    Damping(:,i) = zeta;
end

% Plot everything for w/o Aerodynamic Damping
figure
hold on
plot(Velocities,Frequencies(1,:))
plot(Velocities,Frequencies(2,:))
plot(Velocities,Frequencies(3,:))
plot(Velocities,Frequencies(4,:))
title('Frequency vs Flow Speed (w/o Aerodynamic Damping)')
xlabel('Flow Speed')
ylabel('Frequency')
grid('on')
hold off

figure
hold on
plot(Velocities,100.*Damping(1,:))
plot(Velocities,100.*Damping(2,:))
plot(Velocities,100.*Damping(3,:))
plot(Velocities,100.*Damping(4,:))
title('Damping vs Flow Speed (w/o Aerodynamic Damping)')
xlabel('Flow Speed')
ylabel('Damping (%)')
grid('on')
hold off

% Solve for V_D
V_D = double(solve(det(K)==0,V))
%% Problem No.2 Part b
clear;clc;close all

% Define the variables
s = 7.5;                 % [m]
c = 2;                   % [m]
xf = 0.48*c;             % [m]
m = 200;                 % [kg/m^2]
EI = 2e7;                % [Nm^2]
GJ = 2e6;                % [Nm^2]
aw = 2*pi;               % Lift curve
rho = 1.225;             % [kg/m^3]
ab = xf - c/2;
e = (c/4 + ab)/c;
M_thetad = -1.2;         % Non-dimensional pitch damping derivative

syms V lamda
q = (1/2)*rho*V^2;       % Dynamic pressure

% Set up matrices (Mx_ddot + Dx_dot + Kx = 0)
M = m.*[(s/5)*c   (s/6)*c   (s/4)*((1/2)*c^2 - c*xf)   (s/5)*((1/2)*c^2 - c*xf);
        (s/6)*c   (s/7)*c   (s/5)*((1/2)*c^2 - c*xf)   (s/6)*((1/2)*c^2 - c*xf);
        (s/4)*((1/2)*c^2 - c*xf)   (s/5)*((1/2)*c^2 - c*xf)   (s/3)*((1/3)*c^3 - c^2*xf + c*xf^2)   (s/4)*((1/3)*c^3 - c^2*xf + c*xf^2);
        (s/5)*((1/2)*c^2 - c*xf)   (s/6)*((1/2)*c^2 - c*xf)   (s/4)*((1/3)*c^3 - c^2*xf + c*xf^2)   (s/5)*((1/3)*c^3 - c^2*xf + c*xf^2)];

K = [(4/s^3)*EI   (6/s^3)*EI   (s/4)*q*c*aw   (s/5)*q*c*aw;
     (6/s^3)*EI  (12/s^3)*EI   (s/5)*q*c*aw   (s/6)*q*c*aw;
            0           0       -(s/3)*q*c^2*aw*e + (1/s)*GJ   -(s/4)*q*c^2*aw*e + (1/s)*GJ;
            0           0       -(s/4)*q*c^2*aw*e + (1/s)*GJ   -(s/5)*q*c^2*aw*e + (4/(3*s))*GJ];

D = [(s/(5*V))*(q*c*aw)         (s/(6*V))*(q*c*aw)       0                           0;
     (s/(6*V))*(q*c*aw)         (s/(7*V))*(q*c*aw)       0                           0;
     -(s/(4*V))*(q*c^2*e*aw)   -(s/(5*V))*(q*c^2*e*aw)   -M_thetad*c^3*q*(s/(12*V))  -M_thetad*c^3*q*(s/(16*V));
     -(s/(5*V))*(q*c^2*e*aw)   -(s/(6*V))*(q*c^2*e*aw)   -M_thetad*c^3*q*(s/(16*V))  -M_thetad*c^3*q*(s/(20*V))];

% Reduce order of equations
A = [zeros(4,4) eye(4);
     -inv(M)*K  -inv(M)*D];

% Plug in velocities
Velocities = 0:1:180;
Frequencies = zeros(8,length(Velocities));
Damping = zeros(8,length(Velocities));
for i = 1:length(Velocities)
    v = Velocities(i);
    Temp = double(subs(A,V,v));
    [wn,zeta] = damp(Temp);
    Frequencies(:,i) = wn;
    Damping(:,i) = zeta;
end

% Plot everything for w/Aerodynamic Damping
figure
hold on
plot(Velocities,Frequencies(1,:))
plot(Velocities,Frequencies(2,:))
plot(Velocities,Frequencies(3,:))
plot(Velocities,Frequencies(4,:))
plot(Velocities,Frequencies(5,:))
plot(Velocities,Frequencies(6,:))
plot(Velocities,Frequencies(7,:))
plot(Velocities,Frequencies(8,:))
grid('on')
title('Frequency vs Flow Speed (w/Aerodynamic Damping)')
xlabel('Flow Speed')
ylabel('Frequency')
hold off

figure
hold on
plot(Velocities,100.*Damping(1,:))
plot(Velocities,100.*Damping(2,:))
plot(Velocities,100.*Damping(3,:))
plot(Velocities,100.*Damping(4,:))
plot(Velocities,100.*Damping(5,:))
plot(Velocities,100.*Damping(6,:))
plot(Velocities,100.*Damping(7,:))
plot(Velocities,100.*Damping(8,:))
title('Damping vs Flow Speed (w/Aerodynamic Damping)')
xlabel('Flow Speed')
ylabel('Damping (%)')
grid('on')
hold off

% Reduce order of equations
A = [zeros(4,4) eye(4);
     -inv(M)*K  zeros(4,4)];

% Plug in velocities
Velocities = 0:1:180;
Frequencies = zeros(8,length(Velocities));
Damping = zeros(8,length(Velocities));
for i = 1:length(Velocities)
    v = Velocities(i);
    Temp = double(subs(A,V,v));
    [wn,zeta] = damp(Temp);
    Frequencies(:,i) = wn;
    Damping(:,i) = zeta;
end

% Plot everything for w/o Aerodynamic Damping
figure
hold on
plot(Velocities,Frequencies(1,:))
plot(Velocities,Frequencies(2,:))
plot(Velocities,Frequencies(3,:))
plot(Velocities,Frequencies(4,:))
plot(Velocities,Frequencies(5,:))
plot(Velocities,Frequencies(6,:))
plot(Velocities,Frequencies(7,:))
plot(Velocities,Frequencies(8,:))
title('Frequency vs Flow Speed (w/o Aerodynamic Damping)')
xlabel('Flow Speed')
ylabel('Frequency')
grid('on')
hold off

figure
hold on
plot(Velocities,100.*Damping(1,:))
plot(Velocities,100.*Damping(2,:))
plot(Velocities,100.*Damping(3,:))
plot(Velocities,100.*Damping(4,:))
plot(Velocities,100.*Damping(5,:))
plot(Velocities,100.*Damping(6,:))
plot(Velocities,100.*Damping(7,:))
plot(Velocities,100.*Damping(8,:))
title('Damping vs Flow Speed (w/o Aerodynamic Damping)')
xlabel('Flow Speed')
ylabel('Damping (%)')
grid('on')
hold off

% Solve for V_D
V_D = double(solve(det(K)==0,V))