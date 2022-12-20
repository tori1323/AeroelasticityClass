%% Final Exam - Aeroelasticity - ME597/AAE556
% Victoria Nagorski - 12/12/22

%% Problem No. 1 Part a
clear;clc;close all;
syms x a omega rho EI lambda omega1 omega2 Omega c1

% Solve for the integrals
A1 = int(sin(pi*x/a)^2,x,0,a);
B1 = int(sin(pi*x/a)*sin(2*pi*x/a),x,0,a);
C1 = int(sin(pi*x/a)*cos(pi*x/a),x,0,a);
D1 = int(sin(pi*x/a)*cos(2*pi*x/a),x,0,a);

A2 = int(sin(2*pi*x/a)*sin(pi*x/a),x,0,a);
B2 = int(sin(2*pi*x/a)^2,x,0,a);
C2 = int(sin(2*pi*x/a)*cos(pi*x/a),x,0,a);
D2 = int(sin(2*pi*x/a)*cos(2*pi*x/a),x,0,a);

% Create the matrix
Temp1 = -omega^2*rho/EI + omega1^2;
Temp2 = -omega^2*rho/EI + omega2^2;
Matrix = [(Temp1*A1 + lambda*(pi/a)*C1) (Temp2*B1 + lambda*(2*pi/a)*D1);
          (Temp1*A2 + lambda*(pi/a)*C2) (Temp2*B2 + lambda*(2*pi/a)*D2)]

A = [(a/2)*omega1^2 -(4/3)*lambda;
     (4/3)*lambda    (a/2)*omega2^2];
determ = det(A - eye(2)*Omega);
OMEGA = solve(determ == 0,Omega)
Small_Omega = OMEGA(1);

% Guess
c2 = 1;

% Ax = Omega*x (Solve for eigenvector)
left = A*[c1;c2];
right = Small_Omega*[c1;c2];

% Solve for c1
C1 = solve(left(1) == right(1), c1);

% Check bottom
Left = A*[C1;c2];
Right = Small_Omega*[C1;c2];

Bottom = double(simplify(Left(2)-Right(2)))
if Bottom == 0
    fprintf('c1 and c2 are correct\n')
else
    fprintf('Not correct\n')
end

%% Problem No. 1 Part d
clear;clc;close all;
syms x a omega EI lambda omega1 omega2 omega3 rho N Nx OMEGA Omega real
% Solve for the matrix
n = 3;

% Define equation
start = (-omega^2*rho/EI + (pi*N/a)^4 - Nx/EI*(pi*N/a)^2)*sin(pi*N*x/a) + ...
    (pi*N/a)*lambda*cos(pi*N*x/a);

% Loop over the equation
eqs = sym(zeros(n,n));                  % Initialize the equation matrix
for int_vert = 1:n
    temp = sin(int_vert*pi*x/a);        % Define new sine for new line
    for int_hor = 1:n
        temp2 = subs(start,N,int_hor);
        temp3 = int(temp*temp2,x,0,a);  % Define each c_N part
        eqs(int_vert,int_hor) = temp3;
    end
end

% Clean up A matrix
A = eqs;
A(1,1) =  (a/2)*omega1^2;
A(2,2) =  (a/2)*omega2^2;
A(3,3) =  (a/2)*omega3^2;

% Solve for derminant
determ = det(A-(eye(n).*Omega))
warning('off','all')            % Keep publish clean

% Solve for eigenvalues
OMEGA = solve(determ == 0,Omega,"MaxDegree",3);

% Set two eigenvalues equal (where they collide)
LAMBDA = solve(OMEGA(2)==OMEGA(3),lambda,"MaxDegree",3);
warning('on','all')

% Pull out eigenvalues without imaginary numbers in them
LAMBDA3 = LAMBDA(3);
LAMBDA4 = LAMBDA(4);

% Make sure they don't become imaginary when square-rooted
test1 = sqrt(LAMBDA3);
test2 = sqrt(LAMBDA4);

string(test1);
strfind(ans,'i');

string(test2);
temp = strfind(ans,'i');

% Make sure no imaginary numbers
if isempty(temp) == 1
    fprintf('Velocity is not imaginary\n')
else
    fprintf('Velocity is imaginary\n')
end

fprintf('The lambda term is:\n')
lambda = string((-(75*((11428703*a^3*omega1^6 - 1191016*a^3*omega2^6 + ...
    2703959*a^3*omega3^6 + 4617996*a^3*omega1^2*omega2^4 - ...
    17559642*a^3*omega1^4*omega2^2 + 3784821*a^3*omega1^2*omega3^4 - ...
    16726467*a^3*omega1^4*omega3^2 - 11896698*a^3*omega2^2*omega3^4 - ...
    1044948*a^3*omega2^4*omega3^2 + 28620*3^(1/2)*abs(omega1^2 - ...
    omega3^2)*abs(a)^3*(50463*omega1^8 - 150845*omega1^6*omega2^2 - ...
    51007*omega1^6*omega3^2 + 144273*omega1^4*omega2^4 + ...
    163989*omega1^4*omega2^2*omega3^2 - 5484*omega1^4*omega3^4 - ...
    52152*omega1^2*omega2^6 - 132090*omega1^2*omega2^4*omega3^2 - ...
    31899*omega1^2*omega2^2*omega3^4 + 14289*omega1^2*omega3^6 + ...
    11236*omega2^8 + 7208*omega2^6*omega3^2 + 55233*omega2^4*omega3^4 - ...
    26189*omega2^2*omega3^6 + 2975*omega3^8)^(1/2) + ...
    25883292*a^3*omega1^2*omega2^2*omega3^2)^(4/3) - ...
    2779189261*a^4*omega1^8 - 126247696*a^4*omega2^8 + ...
    168568979*a^4*omega3^8 + 11363*a^2*omega1^4*(11428703*a^3*omega1^6 - ...
    1191016*a^3*omega2^6 + 2703959*a^3*omega3^6 + 4617996*a^3*omega1^2*omega2^4 -...
    17559642*a^3*omega1^4*omega2^2 + 3784821*a^3*omega1^2*omega3^4 - ...
    16726467*a^3*omega1^4*omega3^2 - 11896698*a^3*omega2^2*omega3^4 - ...
    1044948*a^3*omega2^4*omega3^2 + 28620*3^(1/2)*abs(omega1^2 - omega3^2)*abs(a)^3*(50463*omega1^8 - 150845*omega1^6*omega2^2 - 51007*omega1^6*omega3^2 + 144273*omega1^4*omega2^4 + 163989*omega1^4*omega2^2*omega3^2 - 5484*omega1^4*omega3^4 - 52152*omega1^2*omega2^6 - 132090*omega1^2*omega2^4*omega3^2 - 31899*omega1^2*omega2^2*omega3^4 + 14289*omega1^2*omega3^6 + 11236*omega2^8 + 7208*omega2^6*omega3^2 + 55233*omega2^4*omega3^4 - 26189*omega2^2*omega3^6 + 2975*omega3^8)^(1/2) + 25883292*a^3*omega1^2*omega2^2*omega3^2)^(2/3) - 11236*a^2*omega2^4*(11428703*a^3*omega1^6 - 1191016*a^3*omega2^6 + 2703959*a^3*omega3^6 + 4617996*a^3*omega1^2*omega2^4 - 17559642*a^3*omega1^4*omega2^2 + 3784821*a^3*omega1^2*omega3^4 - 16726467*a^3*omega1^4*omega3^2 - 11896698*a^3*omega2^2*omega3^4 - 1044948*a^3*omega2^4*omega3^2 + 28620*3^(1/2)*abs(omega1^2 - omega3^2)*abs(a)^3*(50463*omega1^8 - 150845*omega1^6*omega2^2 - 51007*omega1^6*omega3^2 + 144273*omega1^4*omega2^4 + 163989*omega1^4*omega2^2*omega3^2 - 5484*omega1^4*omega3^4 - 52152*omega1^2*omega2^6 - 132090*omega1^2*omega2^4*omega3^2 - 31899*omega1^2*omega2^2*omega3^4 + 14289*omega1^2*omega3^6 + 11236*omega2^8 + 7208*omega2^6*omega3^2 + 55233*omega2^4*omega3^4 - 26189*omega2^2*omega3^6 + 2975*omega3^8)^(1/2) + 25883292*a^3*omega1^2*omega2^2*omega3^2)^(2/3) - 42061*a^2*omega3^4*(11428703*a^3*omega1^6 - 1191016*a^3*omega2^6 + 2703959*a^3*omega3^6 + 4617996*a^3*omega1^2*omega2^4 - 17559642*a^3*omega1^4*omega2^2 + 3784821*a^3*omega1^2*omega3^4 - 16726467*a^3*omega1^4*omega3^2 - 11896698*a^3*omega2^2*omega3^4 - 1044948*a^3*omega2^4*omega3^2 + 28620*3^(1/2)*abs(omega1^2 - omega3^2)*abs(a)^3*(50463*omega1^8 - 150845*omega1^6*omega2^2 - 51007*omega1^6*omega3^2 + 144273*omega1^4*omega2^4 + 163989*omega1^4*omega2^2*omega3^2 - 5484*omega1^4*omega3^4 - 52152*omega1^2*omega2^6 - 132090*omega1^2*omega2^4*omega3^2 - 31899*omega1^2*omega2^2*omega3^4 + 14289*omega1^2*omega3^6 + 11236*omega2^8 + 7208*omega2^6*omega3^2 + 55233*omega2^4*omega3^4 - 26189*omega2^2*omega3^6 + 2975*omega3^8)^(1/2) + 25883292*a^3*omega1^2*omega2^2*omega3^2)^(2/3) - 5142706*a^3*omega1^6*(11428703*a^3*omega1^6 - 1191016*a^3*omega2^6 + 2703959*a^3*omega3^6 + 4617996*a^3*omega1^2*omega2^4 - 17559642*a^3*omega1^4*omega2^2 + 3784821*a^3*omega1^2*omega3^4 - 16726467*a^3*omega1^4*omega3^2 - 11896698*a^3*omega2^2*omega3^4 - 1044948*a^3*omega2^4*omega3^2 + 28620*3^(1/2)*abs(omega1^2 - omega3^2)*abs(a)^3*(50463*omega1^8 - 150845*omega1^6*omega2^2 - 51007*omega1^6*omega3^2 + 144273*omega1^4*omega2^4 + 163989*omega1^4*omega2^2*omega3^2 - 5484*omega1^4*omega3^4 - 52152*omega1^2*omega2^6 - 132090*omega1^2*omega2^4*omega3^2 - 31899*omega1^2*omega2^2*omega3^4 + 14289*omega1^2*omega3^6 + 11236*omega2^8 + 7208*omega2^6*omega3^2 + 55233*omega2^4*omega3^4 - 26189*omega2^2*omega3^6 + 2975*omega3^8)^(1/2) + 25883292*a^3*omega1^2*omega2^2*omega3^2)^(1/3) + 2382032*a^3*omega2^6*(11428703*a^3*omega1^6 - 1191016*a^3*omega2^6 + 2703959*a^3*omega3^6 + 4617996*a^3*omega1^2*omega2^4 - 17559642*a^3*omega1^4*omega2^2 + 3784821*a^3*omega1^2*omega3^4 - 16726467*a^3*omega1^4*omega3^2 - 11896698*a^3*omega2^2*omega3^4 - 1044948*a^3*omega2^4*omega3^2 + 28620*3^(1/2)*abs(omega1^2 - omega3^2)*abs(a)^3*(50463*omega1^8 - 150845*omega1^6*omega2^2 - 51007*omega1^6*omega3^2 + 144273*omega1^4*omega2^4 + 163989*omega1^4*omega2^2*omega3^2 - 5484*omega1^4*omega3^4 - 52152*omega1^2*omega2^6 - 132090*omega1^2*omega2^4*omega3^2 - 31899*omega1^2*omega2^2*omega3^4 + 14289*omega1^2*omega3^6 + 11236*omega2^8 + 7208*omega2^6*omega3^2 + 55233*omega2^4*omega3^4 - 26189*omega2^2*omega3^6 + 2975*omega3^8)^(1/2) + 25883292*a^3*omega1^2*omega2^2*omega3^2)^(1/3) + 59582*a^3*omega3^6*(11428703*a^3*omega1^6 - 1191016*a^3*omega2^6 + 2703959*a^3*omega3^6 + 4617996*a^3*omega1^2*omega2^4 - 17559642*a^3*omega1^4*omega2^2 + 3784821*a^3*omega1^2*omega3^4 - 16726467*a^3*omega1^4*omega3^2 - 11896698*a^3*omega2^2*omega3^4 - 1044948*a^3*omega2^4*omega3^2 + 28620*3^(1/2)*abs(omega1^2 - omega3^2)*abs(a)^3*(50463*omega1^8 - 150845*omega1^6*omega2^2 - 51007*omega1^6*omega3^2 + 144273*omega1^4*omega2^4 + 163989*omega1^4*omega2^2*omega3^2 - 5484*omega1^4*omega3^4 - 52152*omega1^2*omega2^6 - 132090*omega1^2*omega2^4*omega3^2 - 31899*omega1^2*omega2^2*omega3^4 + 14289*omega1^2*omega3^6 + 11236*omega2^8 + 7208*omega2^6*omega3^2 + 55233*omega2^4*omega3^4 - 26189*omega2^2*omega3^6 + 2975*omega3^8)^(1/2) + 25883292*a^3*omega1^2*omega2^2*omega3^2)^(1/3) + 652676768*a^4*omega1^2*omega2^6 - 3722644104*a^4*omega1^4*omega2^4 + 6143973272*a^4*omega1^6*omega2^2 - 522551332*a^4*omega1^2*omega3^6 - 1965859854*a^4*omega1^4*omega3^4 + 4972783772*a^4*omega1^6*omega3^2 - 151724584*a^4*omega2^2*omega3^6 - 2522099976*a^4*omega2^4*omega3^4 - 147685984*a^4*omega2^6*omega3^2 - 42188*a^2*omega1^2*omega2^2*(11428703*a^3*omega1^6 - 1191016*a^3*omega2^6 + 2703959*a^3*omega3^6 + 4617996*a^3*omega1^2*omega2^4 - 17559642*a^3*omega1^4*omega2^2 + 3784821*a^3*omega1^2*omega3^4 - 16726467*a^3*omega1^4*omega3^2 - 11896698*a^3*omega2^2*omega3^4 - 1044948*a^3*omega2^4*omega3^2 + 28620*3^(1/2)*abs(omega1^2 - omega3^2)*abs(a)^3*(50463*omega1^8 - 150845*omega1^6*omega2^2 - 51007*omega1^6*omega3^2 + 144273*omega1^4*omega2^4 + 163989*omega1^4*omega2^2*omega3^2 - 5484*omega1^4*omega3^4 - 52152*omega1^2*omega2^6 - 132090*omega1^2*omega2^4*omega3^2 - 31899*omega1^2*omega2^2*omega3^4 + 14289*omega1^2*omega3^6 + 11236*omega2^8 + 7208*omega2^6*omega3^2 + 55233*omega2^4*omega3^4 - 26189*omega2^2*omega3^6 + 2975*omega3^8)^(1/2) + 25883292*a^3*omega1^2*omega2^2*omega3^2)^(2/3) + 19462*a^2*omega1^2*omega3^2*(11428703*a^3*omega1^6 - 1191016*a^3*omega2^6 + 2703959*a^3*omega3^6 + 4617996*a^3*omega1^2*omega2^4 - 17559642*a^3*omega1^4*omega2^2 + 3784821*a^3*omega1^2*omega3^4 - 16726467*a^3*omega1^4*omega3^2 - 11896698*a^3*omega2^2*omega3^4 - 1044948*a^3*omega2^4*omega3^2 + 28620*3^(1/2)*abs(omega1^2 - omega3^2)*abs(a)^3*(50463*omega1^8 - 150845*omega1^6*omega2^2 - 51007*omega1^6*omega3^2 + 144273*omega1^4*omega2^4 + 163989*omega1^4*omega2^2*omega3^2 - 5484*omega1^4*omega3^4 - 52152*omega1^2*omega2^6 - 132090*omega1^2*omega2^4*omega3^2 - 31899*omega1^2*omega2^2*omega3^4 + 14289*omega1^2*omega3^6 + 11236*omega2^8 + 7208*omega2^6*omega3^2 + 55233*omega2^4*omega3^4 - 26189*omega2^2*omega3^6 + 2975*omega3^8)^(1/2) + 25883292*a^3*omega1^2*omega2^2*omega3^2)^(2/3) + 64660*a^2*omega2^2*omega3^2*(11428703*a^3*omega1^6 - 1191016*a^3*omega2^6 + 2703959*a^3*omega3^6 + 4617996*a^3*omega1^2*omega2^4 - 17559642*a^3*omega1^4*omega2^2 + 3784821*a^3*omega1^2*omega3^4 - 16726467*a^3*omega1^4*omega3^2 - 11896698*a^3*omega2^2*omega3^4 - 1044948*a^3*omega2^4*omega3^2 + 28620*3^(1/2)*abs(omega1^2 - omega3^2)*abs(a)^3*(50463*omega1^8 - 150845*omega1^6*omega2^2 - 51007*omega1^6*omega3^2 + 144273*omega1^4*omega2^4 + 163989*omega1^4*omega2^2*omega3^2 - 5484*omega1^4*omega3^4 - 52152*omega1^2*omega2^6 - 132090*omega1^2*omega2^4*omega3^2 - 31899*omega1^2*omega2^2*omega3^4 + 14289*omega1^2*omega3^6 + 11236*omega2^8 + 7208*omega2^6*omega3^2 + 55233*omega2^4*omega3^4 - 26189*omega2^2*omega3^6 + 2975*omega3^8)^(1/2) + 25883292*a^3*omega1^2*omega2^2*omega3^2)^(2/3) - 9235992*a^3*omega1^2*omega2^4*(11428703*a^3*omega1^6 - 1191016*a^3*omega2^6 + 2703959*a^3*omega3^6 + 4617996*a^3*omega1^2*omega2^4 - 17559642*a^3*omega1^4*omega2^2 + 3784821*a^3*omega1^2*omega3^4 - 16726467*a^3*omega1^4*omega3^2 - 11896698*a^3*omega2^2*omega3^4 - 1044948*a^3*omega2^4*omega3^2 + 28620*3^(1/2)*abs(omega1^2 - omega3^2)*abs(a)^3*(50463*omega1^8 - 150845*omega1^6*omega2^2 - 51007*omega1^6*omega3^2 + 144273*omega1^4*omega2^4 + 163989*omega1^4*omega2^2*omega3^2 - 5484*omega1^4*omega3^4 - 52152*omega1^2*omega2^6 - 132090*omega1^2*omega2^4*omega3^2 - 31899*omega1^2*omega2^2*omega3^4 + 14289*omega1^2*omega3^6 + 11236*omega2^8 + 7208*omega2^6*omega3^2 + 55233*omega2^4*omega3^4 - 26189*omega2^2*omega3^6 + 2975*omega3^8)^(1/2) + 25883292*a^3*omega1^2*omega2^2*omega3^2)^(1/3) + 11937084*a^3*omega1^4*omega2^2*(11428703*a^3*omega1^6 - 1191016*a^3*omega2^6 + 2703959*a^3*omega3^6 + 4617996*a^3*omega1^2*omega2^4 - 17559642*a^3*omega1^4*omega2^2 + 3784821*a^3*omega1^2*omega3^4 - 16726467*a^3*omega1^4*omega3^2 - 11896698*a^3*omega2^2*omega3^4 - 1044948*a^3*omega2^4*omega3^2 + 28620*3^(1/2)*abs(omega1^2 - omega3^2)*abs(a)^3*(50463*omega1^8 - 150845*omega1^6*omega2^2 - 51007*omega1^6*omega3^2 + 144273*omega1^4*omega2^4 + 163989*omega1^4*omega2^2*omega3^2 - 5484*omega1^4*omega3^4 - 52152*omega1^2*omega2^6 - 132090*omega1^2*omega2^4*omega3^2 - 31899*omega1^2*omega2^2*omega3^4 + 14289*omega1^2*omega3^6 + 11236*omega2^8 + 7208*omega2^6*omega3^2 + 55233*omega2^4*omega3^4 - 26189*omega2^2*omega3^6 + 2975*omega3^8)^(1/2) + 25883292*a^3*omega1^2*omega2^2*omega3^2)^(1/3) - 789942*a^3*omega1^2*omega3^4*(11428703*a^3*omega1^6 - 1191016*a^3*omega2^6 + 2703959*a^3*omega3^6 + 4617996*a^3*omega1^2*omega2^4 - 17559642*a^3*omega1^4*omega2^2 + 3784821*a^3*omega1^2*omega3^4 - 16726467*a^3*omega1^4*omega3^2 - 11896698*a^3*omega2^2*omega3^4 - 1044948*a^3*omega2^4*omega3^2 + 28620*3^(1/2)*abs(omega1^2 - omega3^2)*abs(a)^3*(50463*omega1^8 - 150845*omega1^6*omega2^2 - 51007*omega1^6*omega3^2 + 144273*omega1^4*omega2^4 + 163989*omega1^4*omega2^2*omega3^2 - 5484*omega1^4*omega3^4 - 52152*omega1^2*omega2^6 - 132090*omega1^2*omega2^4*omega3^2 - 31899*omega1^2*omega2^2*omega3^4 + 14289*omega1^2*omega3^6 + 11236*omega2^8 + 7208*omega2^6*omega3^2 + 55233*omega2^4*omega3^4 - 26189*omega2^2*omega3^6 + 2975*omega3^8)^(1/2) + 25883292*a^3*omega1^2*omega2^2*omega3^2)^(1/3) + 3491034*a^3*omega1^4*omega3^2*(11428703*a^3*omega1^6 - 1191016*a^3*omega2^6 + 2703959*a^3*omega3^6 + 4617996*a^3*omega1^2*omega2^4 - 17559642*a^3*omega1^4*omega2^2 + 3784821*a^3*omega1^2*omega3^4 - 16726467*a^3*omega1^4*omega3^2 - 11896698*a^3*omega2^2*omega3^4 - 1044948*a^3*omega2^4*omega3^2 + 28620*3^(1/2)*abs(omega1^2 - omega3^2)*abs(a)^3*(50463*omega1^8 - 150845*omega1^6*omega2^2 - 51007*omega1^6*omega3^2 + 144273*omega1^4*omega2^4 + 163989*omega1^4*omega2^2*omega3^2 - 5484*omega1^4*omega3^4 - 52152*omega1^2*omega2^6 - 132090*omega1^2*omega2^4*omega3^2 - 31899*omega1^2*omega2^2*omega3^4 + 14289*omega1^2*omega3^6 + 11236*omega2^8 + 7208*omega2^6*omega3^2 + 55233*omega2^4*omega3^4 - 26189*omega2^2*omega3^6 + 2975*omega3^8)^(1/2) + 25883292*a^3*omega1^2*omega2^2*omega3^2)^(1/3) + 611196*a^3*omega2^2*omega3^4*(11428703*a^3*omega1^6 - 1191016*a^3*omega2^6 + 2703959*a^3*omega3^6 + 4617996*a^3*omega1^2*omega2^4 - 17559642*a^3*omega1^4*omega2^2 + 3784821*a^3*omega1^2*omega3^4 - 16726467*a^3*omega1^4*omega3^2 - 11896698*a^3*omega2^2*omega3^4 - 1044948*a^3*omega2^4*omega3^2 + 28620*3^(1/2)*abs(omega1^2 - omega3^2)*abs(a)^3*(50463*omega1^8 - 150845*omega1^6*omega2^2 - 51007*omega1^6*omega3^2 + 144273*omega1^4*omega2^4 + 163989*omega1^4*omega2^2*omega3^2 - 5484*omega1^4*omega3^4 - 52152*omega1^2*omega2^6 - 132090*omega1^2*omega2^4*omega3^2 - 31899*omega1^2*omega2^2*omega3^4 + 14289*omega1^2*omega3^6 + 11236*omega2^8 + 7208*omega2^6*omega3^2 + 55233*omega2^4*omega3^4 - 26189*omega2^2*omega3^6 + 2975*omega3^8)^(1/2) + 25883292*a^3*omega1^2*omega2^2*omega3^2)^(1/3) + 2089896*a^3*omega2^4*omega3^2*(11428703*a^3*omega1^6 - 1191016*a^3*omega2^6 + 2703959*a^3*omega3^6 + 4617996*a^3*omega1^2*omega2^4 - 17559642*a^3*omega1^4*omega2^2 + 3784821*a^3*omega1^2*omega3^4 - 16726467*a^3*omega1^4*omega3^2 - 11896698*a^3*omega2^2*omega3^4 - 1044948*a^3*omega2^4*omega3^2 + 28620*3^(1/2)*abs(omega1^2 - omega3^2)*abs(a)^3*(50463*omega1^8 - 150845*omega1^6*omega2^2 - 51007*omega1^6*omega3^2 + 144273*omega1^4*omega2^4 + 163989*omega1^4*omega2^2*omega3^2 - 5484*omega1^4*omega3^4 - 52152*omega1^2*omega2^6 - 132090*omega1^2*omega2^4*omega3^2 - 31899*omega1^2*omega2^2*omega3^4 + 14289*omega1^2*omega3^6 + 11236*omega2^8 + 7208*omega2^6*omega3^2 + 55233*omega2^4*omega3^4 - 26189*omega2^2*omega3^6 + 2975*omega3^8)^(1/2) + 25883292*a^3*omega1^2*omega2^2*omega3^2)^(1/3) + 5499373704*a^4*omega1^2*omega2^2*omega3^4 + 5487257904*a^4*omega1^2*omega2^4*omega3^2 - 10986631608*a^4*omega1^4*omega2^2*omega3^2 - 5402184*a^3*omega1^2*omega2^2*omega3^2*(11428703*a^3*omega1^6 - 1191016*a^3*omega2^6 + 2703959*a^3*omega3^6 + 4617996*a^3*omega1^2*omega2^4 - 17559642*a^3*omega1^4*omega2^2 + 3784821*a^3*omega1^2*omega3^4 - 16726467*a^3*omega1^4*omega3^2 - 11896698*a^3*omega2^2*omega3^4 - 1044948*a^3*omega2^4*omega3^2 + 28620*3^(1/2)*abs(omega1^2 - omega3^2)*abs(a)^3*(50463*omega1^8 - 150845*omega1^6*omega2^2 - 51007*omega1^6*omega3^2 + 144273*omega1^4*omega2^4 + 163989*omega1^4*omega2^2*omega3^2 - 5484*omega1^4*omega3^4 - 52152*omega1^2*omega2^6 - 132090*omega1^2*omega2^4*omega3^2 - 31899*omega1^2*omega2^2*omega3^4 + 14289*omega1^2*omega3^6 + 11236*omega2^8 + 7208*omega2^6*omega3^2 + 55233*omega2^4*omega3^4 - 26189*omega2^2*omega3^6 + 2975*omega3^8)^(1/2) + 25883292*a^3*omega1^2*omega2^2*omega3^2)^(1/3) - 7841880*3^(1/2)*a*omega1^2*abs(omega1^2 - omega3^2)*abs(a)^3*(50463*omega1^8 - 150845*omega1^6*omega2^2 - 51007*omega1^6*omega3^2 + 144273*omega1^4*omega2^4 + 163989*omega1^4*omega2^2*omega3^2 - 5484*omega1^4*omega3^4 - 52152*omega1^2*omega2^6 - 132090*omega1^2*omega2^4*omega3^2 - 31899*omega1^2*omega2^2*omega3^4 + 14289*omega1^2*omega3^6 + 11236*omega2^8 + 7208*omega2^6*omega3^2 + 55233*omega2^4*omega3^4 - 26189*omega2^2*omega3^6 + 2975*omega3^8)^(1/2) + 6067440*3^(1/2)*a*omega2^2*abs(omega1^2 - omega3^2)*abs(a)^3*(50463*omega1^8 - 150845*omega1^6*omega2^2 - 51007*omega1^6*omega3^2 + 144273*omega1^4*omega2^4 + 163989*omega1^4*omega2^2*omega3^2 - 5484*omega1^4*omega3^4 - 52152*omega1^2*omega2^6 - 132090*omega1^2*omega2^4*omega3^2 - 31899*omega1^2*omega2^2*omega3^4 + 14289*omega1^2*omega3^6 + 11236*omega2^8 + 7208*omega2^6*omega3^2 + 55233*omega2^4*omega3^4 - 26189*omega2^2*omega3^6 + 2975*omega3^8)^(1/2) + 1774440*3^(1/2)*a*omega3^2*abs(omega1^2 - omega3^2)*abs(a)^3*(50463*omega1^8 - 150845*omega1^6*omega2^2 - 51007*omega1^6*omega3^2 + 144273*omega1^4*omega2^4 + 163989*omega1^4*omega2^2*omega3^2 - 5484*omega1^4*omega3^4 - 52152*omega1^2*omega2^6 - 132090*omega1^2*omega2^4*omega3^2 - 31899*omega1^2*omega2^2*omega3^4 + 14289*omega1^2*omega3^6 + 11236*omega2^8 + 7208*omega2^6*omega3^2 + 55233*omega2^4*omega3^4 - 26189*omega2^2*omega3^6 + 2975*omega3^8)^(1/2)))/(76225024*(11428703*a^3*omega1^6 - 1191016*a^3*omega2^6 + 2703959*a^3*omega3^6 + 4617996*a^3*omega1^2*omega2^4 - 17559642*a^3*omega1^4*omega2^2 + 3784821*a^3*omega1^2*omega3^4 - 16726467*a^3*omega1^4*omega3^2 - 11896698*a^3*omega2^2*omega3^4 - 1044948*a^3*omega2^4*omega3^2 + 28620*3^(1/2)*abs(omega1^2 - omega3^2)*abs(a)^3*(50463*omega1^8 - 150845*omega1^6*omega2^2 - 51007*omega1^6*omega3^2 + 144273*omega1^4*omega2^4 + 163989*omega1^4*omega2^2*omega3^2 - 5484*omega1^4*omega3^4 - 52152*omega1^2*omega2^6 - 132090*omega1^2*omega2^4*omega3^2 - 31899*omega1^2*omega2^2*omega3^4 + 14289*omega1^2*omega3^6 + 11236*omega2^8 + 7208*omega2^6*omega3^2 + 55233*omega2^4*omega3^4 - 26189*omega2^2*omega3^6 + 2975*omega3^8)^(1/2) + 25883292*a^3*omega1^2*omega2^2*omega3^2)^(2/3)))^(1/2)/2)
%% Problem No. 2 Part a
clear;clc;close all;
syms y s h(t) x xf theta(t) t EI GJ aw rho e M_thetad c m dh dtheta p(t) dp M

% Define the variables... The function explains variables better
mult = 1;           
fuse_mass = mult*200*7.5*2;        % M = mult*m_wing
vars = [s;c;xf;m;EI;GJ;aw;rho;e;M_thetad;M];
nums = [7.5;2;0.48*2;200;2e7;2e6;2*pi;1.225;0.23;-1.2;fuse_mass];
bend_mode = (y/s)^2*h;
dB = (y/s)^2*dh;
tor_mode = (x-xf)*(y/s)*theta;
dT = (y/s)*dtheta;
plunge_mode = p;
modes = [h theta p];
dP = dp;
d = [dh dtheta dp];
damping = true;                    % Not set C = 0 with analysis?

[sys1,~] = FlexWing_Lagrange(vars,nums,bend_mode,tor_mode,damping,modes,dB,dT,d,plunge_mode,dP,mult);

% For publishing purposes
M1 = simplify(sys1.M);
K1 = simplify(sys1.K);
C1 = simplify(sys1.C);

display(M1);
display(K1);
display(C1);

% For 10*m_wing
% Define the variables
mult = 10;          % M = mult*m_wing
fuse_mass = mult*200*7.5*2;
nums = [7.5;2;0.48*2;200;2e7;2e6;2*pi;1.225;0.23;-1.2;fuse_mass];
[~,~] = FlexWing_Lagrange(vars,nums,bend_mode,tor_mode,damping,modes,dB,dT,d,plunge_mode,dP,mult);

%% Problem No. 2 Part b
clear;clc;close all;
syms y s h1(t) h2(t) x xf theta1(t) theta2(t) t EI GJ aw rho e...
    M_thetad c m dh1 dh2 dtheta1 dtheta2 p(t) dp M

% Define the variables... The function explains variables better
mult = 1;           % M = mult*m_wing
fuse_mass = mult*200*7.5*2;
vars = [s;c;xf;m;EI;GJ;aw;rho;e;M_thetad;M];
nums = [7.5;2;0.48*2;200;2e7;2e6;2*pi;1.225;0.23;-1.2;fuse_mass];
bend_mode = (y/s)^2*h1 + (y/s)^3*h2;
dB = (y/s)^2*dh1 + (y/s)^3*dh2;
tor_mode = (x-xf)*(y/s)*theta1 + (x-xf)*(y/s)^2*theta2;
dT = (y/s)*dtheta1 + (y/s)^2*dtheta2;
plunge_mode = p;
modes = [h1 h2 theta1 theta2 p];
dP = dp;
d = [dh1 dh2 dtheta1 dtheta2 dp];
damping = true;                    % Not set C = 0 with analysis?

[sys2,~] = FlexWing_Lagrange(vars,nums,bend_mode,tor_mode,damping,modes,dB,dT,d,plunge_mode,dP,mult);

% For publishing purposes
M2 = simplify(sys2.M);
K2 = simplify(sys2.K);
C2 = simplify(sys2.C);

display(M2);
display(K2);
display(C2);

% For 10*m_wing
% Define the variables
mult = 10;          % M = mult*m_wing
fuse_mass = mult*200*7.5*2;
nums = [7.5;2;0.48*2;200;2e7;2e6;2*pi;1.225;0.23;-1.2;fuse_mass];

[~,~] = FlexWing_Lagrange(vars,nums,bend_mode,tor_mode,damping,modes,dB,dT,d,plunge_mode,dP,mult);

%% Dummy section because of local function
clear;clc;close all;

%%  ======================= Local Function ============================
function [sys,other] = FlexWing_Lagrange(var,nums,bend,tor,damping,mode,dB,dT,d,plunge,dP,mult)
% This function is for a flexible wing with damping
% Outputs:
%   - sys: Mass [M], Damper [C], Spring [K], and Plant [A]
%   - other: Damping and Frequencies
% Inputs:
%   - var: variables
%   - nums: values for the variables
%   - bend: bending modes definition
%   - tor: torsional modes definition
%   - damping: True or False whether aerodynamic damping is to be
%              considered
%   - mode: Modes to consider
%   - dB: delta bend
%   - dT: delta theta
%   - d: delta modes
%   - plunge: plunge mode defintion
%   - dP: delta plunge
%
% Start code below::
    syms t m x c y s EI GJ xf V rho q aw M_thetad e M
    % ============== Solve for Lagrange's equations ==============
    % Solve for kinetic energy
    part1 = diff(bend,t);     % Derivative w.r.t time of bending modes
    part2 = diff(tor,t);      % Derivative w.r.t time of torisional modes
    part3 = diff(plunge,t);   % Derivative w.r.t time of plunging modes
    together = (part1+part2+part3)^2;
    FuselageT = (1/2)*M*part3^2;
    T = simplify(m/2*int(int(together,x,0,c),y,0,s) + FuselageT);

    % Solve for potential energy
    U_pt1 = (1/2)*int(diff(bend,y,2)^2,y,0,s)*EI;
    U_pt2 = (1/2)*int(diff(tor/(x-xf),y)^2,y,0,s)*GJ;
    U = simplify(U_pt1 + U_pt2);

    % Solve for delta work
    wPart1 = q*c*aw * (tor/(x-xf) + diff(bend,t)/V + diff(plunge,t)/V)*(dP + dB);
    wPart2 = q*c^2*(e*aw*(tor/(x-xf) + diff(bend,t)/V + diff(plunge,t)/V) + M_thetad*(diff(tor/(x-xf),t)*c/(4*V)))*dT;
    dW = -int(wPart1,y,0,s) + int(wPart2,y,0,s);

    % Solve Lagrange's equation
    mode = mode(t);
    eq = sym(zeros(length(mode),1));
    for i = 1:length(mode)
        Mode = mode(i);
        dd = d(i);
        dot = diff(Mode,t);
        eq(i,1) = simplify(diff(diff(T,dot),t) - diff(T,Mode) + diff(U,Mode) - diff(dW,dd));
    end
    
    % Create the matrices
    sys.M = simplify(jacobian(eq,diff(mode,t,2)));      % Mass
    sys.K = simplify(jacobian(eq,mode));                % Spring
    if damping == false
        sys.C = zeros(length(mode),length(mode));
    else
        sys.C = simplify(jacobian(eq,diff(mode,t)));    % Damper
    end
    sys.A = [zeros(length(mode),length(mode))    eye(length(mode));
                    -inv(sys.M)*sys.K           -inv(sys.M)*sys.C];

    % ============= Solve for Frequency and Damping ==============
    A = simplify(subs(sys.A,q,(1/2)*rho*V^2));
    A = subs(A,var,nums);

    % Plug in velocities
    Velocities = 1:0.1:100;
    other.frequencies = zeros(length(mode)*2,length(Velocities));
    other.damping = zeros(length(mode)*2,length(Velocities));
    for i = 1:length(Velocities)
        v = Velocities(i);
        Temp = double(subs(A,V,v));
        p = eig(Temp);
        wrad = abs(p);
        zeta = -real(p) ./ wrad;
        if imag(p) == 0
            wrad = 0;
        end
        whz = wrad ./ (2*pi);
        
        [whz, orderf] = sort(whz);

        other.frequencies(:,i) = whz;
        other.damping(:,i) = zeta(orderf);
    end
    
    % Plot everything
    figure
    hold on
    for i = 1:(length(mode)*2)
         plot(Velocities,other.frequencies(i,:))
    end
    title(['Frequency vs Flow Speed || M = ',num2str(mult),' kg'])
    xlabel('Flow Speed')
    ylabel('Frequency')
    xlim([0 100])
    ylim([0 6])
    grid('on')
    hold off

    figure
    hold on
    for i = 1:(length(mode)*2)
         plot(Velocities,other.damping(i,:).*100)
    end
    title(['Damping Ratio vs Flow Speed || M = ',num2str(mult),' kg'])
    xlabel('Flow Speed')
    ylabel('Damping Ratio')
    xlim([0 100])
    ylim([-5 10])
    line([0 100],[0 0],'Color','red','LineStyle','--')
    grid('on')
    hold off
end