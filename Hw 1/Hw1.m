%% Hw #1 - Aeroelasticity - ME597/AAE556
% Victoria Nagorski - 9/2/22

%% Problem No. 1c
clear;clc;close all;
% Initialize Variables
theta = 0:0.001:0.5;
K_e = @(qHat) (1-qHat) + 8*qHat.*theta.^2;
q_hat = [0,0.5,1,1.5];

% Graphing
figure
hold on
for i = 1:length(q_hat)
    qHat = q_hat(1,i);          % Pull out the q_hat value being used
    temp_Ke = K_e(qHat);        % Solve for the K_e values
    plot(theta,temp_Ke)         % Graph everything
end
title('$K\_{e}$ vs $\theta$','Interpreter','latex')
xlabel('$\theta$','Interpreter','latex')
ylabel('$K\_{e}$','Interpreter','latex')
legend('q_{hat} = 0','q_{hat} = 0.5','q_{hat} = 1','q_{hat} = 1.5')
hold off

%% Problem No. 1d
clear;clc;close all;

% Initialize Variables
theta = -0.5:0.0001:0.5;
P_hat = @(qHat) (1-qHat).*theta + 8/3*qHat.*theta.^3;
q_hat = [0,0.5,1,1.5];

% Graphing
figure
hold on
for i = 1:length(q_hat)
    qHat = q_hat(1,i);            % Pull out the q_hat value being used
    temp_Phat = P_hat(qHat);      % Solve for the P_hat values
    plot(theta,temp_Phat)         % Graph everything
end
x = [-0.5 0.5];
y = [0 0];
line(x,y,'Color','red','LineStyle','--')
title('$\hat{P}$ vs $\theta$ ','Interpreter','latex')
xlabel('$\theta$[rads]','Interpreter','latex')
ylabel('$\hat{P}$','Interpreter','latex')
legend('q_{hat} = 0','q_{hat} = 0.5','q_{hat} = 1','q_{hat} = 1.5')
hold off

%% Problem No. 2a
clear;clc;close all;
% Initialize Variables
V = 120;            % Velocity [m/s]
rho = 1.112;        % Density [kg/m^3]
W = 70000;          % Total weight [N]
muRatio = 1/2;      % Mu/e
q = (1/2)*rho*V^2;  % Dynamic pressure (typical) [N/m^2]
S = 125;            % Wing  area [m^2]
a1 = 6;             % Alpha vs lift coefficient slope
n = 0:0.0001:1;

% Define the Equations
alpha_0 = @(qHat) (n.*W./(2*q*S*a1)) .* (1-qHat*(1 + muRatio));
q_hat = [0,0.5,1,1.5];

% Graphing
figure
hold on
for i = 1:length(q_hat)
    qHat = q_hat(1,i);              % Pull out the q_hat value being used
    temp_alpha0 = alpha_0(qHat);    % Solve for the P_hat values
    plot(n,temp_alpha0)             % Graph everything
end
title('$\alpha_0$ vs n','Interpreter','latex')
xlabel('n','Interpreter','latex')
ylabel('$\alpha_0$','Interpreter','latex')
legend('q_{hat} = 0','q_{hat} = 0.5','q_{hat} = 1','q_{hat} = 1.5')
hold off

%% Problem No. 2c
clear;clc;close all;
% Initialize Variables
V = 120;            % Velocity [m/s]
rho = 1.112;        % Density [kg/m^3]
W = 70000;          % Total weight [N]
muRatio = 1/2;      % Mu/e
q = (1/2)*rho*V^2;  % Dynamic pressure (typical) [N/m^2]
S = 125;            % Wing  area [m^2]
a1 = 6;             % Alpha vs lift coefficient slope
n = 0:0.0001:1;

% Define the Equations
theta = @(qHat) (n.*W./(2*q*S*a1)) .* (qHat*(1 + muRatio));
q_hat = [0,0.5,1,1.5];

% Graphing
figure
hold on
for i = 1:length(q_hat)
    qHat = q_hat(1,i);             % Pull out the q_hat value being used
    temp_theta = theta(qHat);      % Solve for the P_hat values
    plot(n,temp_theta)             % Graph everything
end
title('$\theta$ vs n','Interpreter','latex')
xlabel('n','Interpreter','latex')
ylabel('$\theta$[rad]','Interpreter','latex')
legend('q_{hat} = 0','q_{hat} = 0.5','q_{hat} = 1','q_{hat} = 1.5')
hold off
