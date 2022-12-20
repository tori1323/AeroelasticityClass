%% Hw #2 - Aeroelasticity - ME597/AAE556
% Victoria Nagorski - 9/19/22

%% Problem 1b
syms q R c M s W y rho

% Solve for lift
P = [4/3*q*R    1/2*q*R 4/15*q*R;
        1/3*R*c*q  (4/35*q*c*R-4/3*M/s) (1/12*q*c*R-M/s);
        4/35*q*c*R (1/12*q*c*R-M/s) (16/315*q*c*R-16/15*M/s)];
thetas = P^-1 * [W;0;0]

dL_dy = simplify(q * R/s * (1 - y^2/s^2) * [1 y/s y^2/s^2]*thetas)


% Solve for V_d
alpha0 = thetas(1);
[~,D] = numden(alpha0)
q = simplify(solve(D,q))

q_A = q(2)
V_A = sqrt(2/rho * q_A)
