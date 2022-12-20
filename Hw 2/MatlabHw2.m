%% Hw #2 - Aeroelasticity - ME597/AAE556
% Victoria Nagorski - 9/19/22

%% Problem No. 3b
clear;clc;close all;
% Initialize Variables
Sr = 0:0.001:1;
q_hat = @(K_r) (((K_r./Sr)+1+K_r)-sqrt(((K_r./Sr)+1+K_r).^2-(4.*K_r./Sr)))./2;
Kr = [1,10];

% Graphing
figure
hold on
for i = 1:length(Kr)
    K_r = Kr(1,i);          % Pull out the q_hat value being used
    temp_q = q_hat(K_r);        % Solve for the K_e values
    plot(Sr,temp_q)         % Graph everything
end
title('$\hat{q}$ vs $S_R$','Interpreter','latex')
xlabel('$S_R$','Interpreter','latex')
ylabel('$\hat{q}$','Interpreter','latex')
legend('k_{R} = 1','k_{R} = 10')
hold off