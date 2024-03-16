% This code is not running in the name 21065100.m if we run in Z.m
% (Z can be any word) then it will be run perfectly


format short
disp("==============================================")
disp('     NAME - SHRIDHAR KUMAR    ')
disp('     ROLL NO.- 21065100     ')
disp('     CE336 ASSIGNEMENT-2 ( MATLAB CODE )     ')
disp("==============================================") 

Elements = input('The No. of elements in the beam = '); %take 3 elements
Nodes = Elements+1;
K = zeros(((Nodes)*2));
E = 5000 * sqrt(25) * 10^3;% E-->modulus of elasticity and 25 is the Grade of concrete
b=.150;
d=0.3;
I=(b*d^3)/12;%moment of inertia
L = input('The length of the span of the beam in colm matrix = ')  %jitne  spans (Nodes-1) hai utne colms or rows ....take L=[3, 3, 3];


for i=1:Elements

     k1=(E*I)*[12/(L(i))^3 6/(L(i))^2 -12/(L(i))^2 6/(L(i))^2; 
        6/(L(i))^2 4/(L(i)) -6/(L(i))^2 2/(L(i)) ;
        -12/(L(i))^3 -6/(L(i))^2 12/(L(i))^3 -6/(L(i))^2;
        6/(L(i))^2 2/(L(i)) -6/(L(i))^2 4/(L(i))];
   K(2*i-1:2*i+2,2*i-1:2*i+2)=K(2*i-1:2*i+2,2*i-1:2*i+2)+k1;
end
fprintf("________The Globle stiffness matrix of the beam______")
  disp(K);
% Calculate the total load(external+internal) vector At
% Applied load P
P=100;
%Load vector Px
Px = [2 *P, P, 4 * P];
%joint load Aj
Aj = [0; 0; 0; P* 3; P; 0; 0; 0];


Ar = zeros(2*Nodes, 1);
for i = 1:Elements
   % Ax--->member load
    Ax1 = Px(i);
    Ax = [Ax1 / 2; Ax1 * L(i) / 8; Ax1  / 2; -Ax1  * L(i) / 8];
   
    Ar(2 * i - 1:2 * i +2, 1) =  Ar(2 * i - 1:2 * i +2, 1) - Ax; % Ar--> equvalent load
end
% Total load vector At
 At = Ar +Aj ;

 % Extract the Aft vector
 % in the beam there is 3 free degree of freedom (X1)
 X1 = 3;
 % which have the values
 X = [4,6,8];  % according to numbering
  % in the beam there is 5 restrained degree of freedom (X1)
 Y1 = 5;
 % which have the values
Y = [1, 2, 3, 5, 7];

Aft = zeros(X1, 1);
for i = 1:X1
    Aft(i, 1) = At(X(i));
end

% for calculating displacement we also need to extract K (stiffness matrix)
Kfd = zeros(X1, X1);
for i = 1:X1
    for j = 1:X1
        Kfd(i, j) = K(X(i), X(j));
    end
end
 % Calculate Df (Displacements of Free DOFs Matrix)
Df = inv(Kfd) * Aft;

fprintf('Df (Displacements of Free DOFs Matrix):\n');
disp(Df);

% Calculate the restrained DOFs and support reactions
Krf = zeros(Y1, X1);
for i = 1:Y1
    for j = 1:X1
        Krf(i, j) = K(Y(i), X(j));
    end
end
% Rse -->restrain support reactions
Rse = zeros(Y1, 1);
for i = 1:Y1
   Rse(i, 1) =At(Y(i));
end
% Rs-->support reactions
Rs = -Rse + Krf * Df;

% Display support reactions
fprintf('Support Reactions :\n');
disp(Rs);

% Calculate member action forces
for mi = 1:Elements 
    Ax = Px(mi);
    Amfi = [Ax / 2; Ax * L(mi) / 8; Ax / 2; -Ax *L(mi) / 8];
    Dmi = zeros(4, 1);

    for i = 1:X1
        if X(i) >= 2 * mi - 1 && X(i) <= 2 * mi + 2
            Dmi(X(i) - (mi - 1) * 2) = Df(i);
        end
    end

    Kmi = (E * I) / L(mi)^3 * ...
        [12, 6 * (mi), -12, 6 * L(mi); ...
         6 * L(mi), 4 * L(mi)^2, -6 * L(mi), 2 * L(mi)^2; ...
         -12, -6 * L(mi), 12, -6 * L(mi); ...
         6 * L(mi), 2 * L(mi)^2, -6 * L(mi), 4 * L(mi)^2];

    Ami = Amfi + Kmi * Dmi;

    % Display Ami as member end reactions
    fprintf('Am%d (Member %d End Reactions):\n', mi, mi);
   disp(Ami);
end








  
