function [sys,vel,other] = FlexWing_Lagrange(var,nums,bend,tor,damping,mode,dB,dT,d)
% This function is for a flexible wing with damping
% Outputs:
%   - sys:
%   - vel:
%   - other:
% Inputs:
%   - var: variables
%   - nums: values for the variables
%   - bend: bending modes definition
%   - tor: torsional modes definition
%   - damping: True or False whether aerodynamic damping is to be
%              considered
%   - mode: Modes to consider
%   - dB: 
%   - dT:
%   - d:
% Start code below::
    syms t m x c y s EI GJ xf V rho q aw M_thetad e
    % ============== Solve for Lagrange's equations ==============
    % Solve for kinetic energy
    part1 = diff(bend,t);     % Derivative w.r.t time of bending modes
    part2 = diff(tor,t);      % Derivative w.r.t time of torisional modes
    together = (part1+part2)^2;
    T = simplify(m/2*int(int(together,x,0,c),y,0,s));

    % Solve for potential energy
    U_pt1 = (1/2)*int(diff(bend,y,2)^2,y,0,s)*EI;
    U_pt2 = (1/2)*int(diff(tor/(x-xf),y)^2,y,0,s)*GJ;
    U = simplify(U_pt1 + U_pt2);

    % Solve for delta work
    wPart1 = q*c*aw * (tor/(x-xf) + diff(bend,t)/V)*dB;
    wPart2 = q*c^2*(e*aw*(tor/(x-xf) + diff(bend,t)/V) + M_thetad*(diff(tor/(x-xf),t)*c/(4*V)))*dT;
    dW = -int(wPart1,y,0,s) + int(wPart2,y,0,s);

    % Solve Lagrange's equation
    mode = mode(t);
    eq = sym(zeros(length(mode),1));
    for i = 1:length(mode)
        M = mode(i);
        dd = d(i);
        dot = diff(M,t);
        eq(i,1) = simplify(diff(diff(T,dot),t) - diff(T,M) + diff(U,M) - diff(dW,dd));
    end
    
    % Create the matrices
    sys.M = jacobian(eq,diff(mode,t,2)); % Mass
    sys.K = jacobian(eq,mode); % 
    if damping == false
        sys.C = zeros(length(mode),length(mode));
    else
        sys.C = jacobian(eq,diff(mode,t));
    end
    sys.A = [zeros(length(mode),length(mode))    eye(length(mode));
                    -inv(sys.M)*sys.K           -inv(sys.M)*sys.C];

    % ============= Solve for Frequency and Damping ==============
    A = simplify(subs(sys.A,q,(1/2)*rho*V^2));
    A = subs(A,var,nums);

    % Plug in velocities
    Velocities = 1:0.1:180;
    other.frequencies = zeros(length(mode)*2,length(Velocities));
    other.damping = zeros(length(mode)*2,length(Velocities));
    for i = 1:length(Velocities)
        v = Velocities(i);
        Temp = double(subs(A,V,v));
        [wn,zeta] = damp(Temp);
        other.frequencies(:,i) = wn./(2*pi);
        other.damping(:,i) = zeta;
%         p = eig(Temp);
%         wrad = abs(p);
%         zeta = -real(p) ./ wrad;
%         if imag(p) == 0
%             wrad = 0;
%         end
%         whz = wrad ./ (2*pi);
%         
%         [whz, orderf] = sort(whz);
% 
%         other.frequencies(:,i) = whz;
%         other.damping(:,i) = zeta(orderf);

    end
    
    

    % Plot everything
    figure
    hold on
    for i = 1:(length(mode)*2)
         plot(Velocities,other.frequencies(i,:))
    end
    title('Frequency vs Flow Speed (w/Aerodynamic Damping)')
    xlabel('Flow Speed')
    ylabel('Frequency')
    grid('on')
    xlim([10 180])
    ylim([2.5 4.5])
    hold off

    figure
    hold on
    for i = 1:(length(mode)*2)
         plot(Velocities,other.damping(i,:).*100)
    end
    title('Damping Ratio vs Flow Speed')
    xlabel('Flow Speed')
    ylabel('Damping Ratio')
    grid('on')
    xlim([10 180])
    ylim([-5 10])
    hold off


    % ================= Solve for the Velocities =================
    % Find where damping crosses 0
    vel.flutter = 0;
    % det(K) = 0
    vel.divergence = double(solve(det(sys.K)==0,V));

end