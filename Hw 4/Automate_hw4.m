syms y s h(t) x xf theta(t) t EI GJ aw rho e M_thetad c m dh dtheta

% Define the variables
vars = [s;c;xf;m;EI;GJ;aw;rho;e;M_thetad];
nums = [7.5;2;0.48*2;200;2e7;2e6;2*pi;1.225;0.2300;-1.2];
% nums = [7.5;2;0.48*2;200;2e7;2e6;2*pi;1.225;0.2300;0];
bend_mode = (y/s)^2*h;
dB = (y/s)^2*dh;
tor_mode = (x-xf)*(y/s)*theta;
dT = (y/s)*dtheta;
modes = [h theta];
d = [dh dtheta];

[sys,vel,other] = FlexWing_Lagrange(vars,nums,bend_mode,tor_mode,false,modes,dB,dT,d);
%%
syms y s h(t) x xf theta(t) t EI GJ aw rho e M_thetad c m dh dtheta

% Define the variables
vars = [s;c;xf;m;EI;GJ;aw;rho;e;M_thetad];
nums = [7.5;2;0.48*2;200;2e7;2e6;2*pi;1.225;0.2300;-1.2];

syms h1(t) h2(t) dh1 dh2 theta1(t) theta2(t) dtheta1 dtheta2
bend_mode = (y/s)^2*h1 + (y/s)^3*h2;
dB = (y/s)^2*dh1 + + (y/s)^3*dh2;
tor_mode = (x-xf)*(y/s)*theta1 + (x-xf)*(y/s)^2*theta2;
dT = (y/s)*dtheta1 + (y/s)^2*dtheta2;
modes = [h1 h2 theta1 theta2];
d = [dh1 dh2 dtheta1 dtheta2];

[sys,vel,other] = FlexWing_Lagrange(vars,nums,bend_mode,tor_mode,true,modes,dB,dT,d);