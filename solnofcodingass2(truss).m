%%*************************************************************************
%*********[2D Frame Analysis by Direct Stiffness Matrix Method]************
%**************************************************************************
clc
clear
close all
%%
scale_factor = 1; % scale factor for deflection 
scale___factor = 1000000; % scale factor for S.F.D 
scale_____factor= 1000000; % scale factor for B.M.D
%% Problem Description 
m = 3; %No of members 
Nj = 4;  %No of joints 
FF = [4 6 7 9 10 12]'; %Free DOFs
BC = [1 2 3 5 8 11]'; %Restrained DOFs
P = 100000; %Value of point load in "N" 
Xj=3;  %Defining the length "m"
jC = [0 0;Xj 0;2*Xj 0;3*Xj 0]; %Joint Coordinates 
A=0.15*0.3*[1,1,1]; %Material Section Properties (area)"m2"
E=(2.5E+10)*[1,1,1]; I = (0.15*0.3*0.3*0.3)*[1/12,1/12,1/12]; %Material Section Properties (Young's Modulus)"N/m2 & area moment of Inertia"
P_m =[-2*P 0 0 0;-P 0 0 0;-4*P 0 0 0]; %Defining load variable P  
m_Nj = [1 2;2 3;3 4]; %Member Starting and end nodes
%% joint stiffness and load vector 
Kj=zeros(3*Nj,3*Nj); %Predefining Joint Stiffness matrix 
A_j = [0 0 0 0 0 P*Xj 0 P 0 0 0 0]'; %Joint load vector 
A_e = zeros(3*Nj,1); %Predefining Equivalent joint load  
A_f = zeros(6,m); %Predefining Fixed-end force matrix
j_d = zeros(2*Nj,1); % Displaced Joint Coordinates   
for c=1:m 
    j=m_Nj(c,1) ; k=m_Nj(c,2) ;
    L = sqrt(((jC(k,1) - jC(j,1))^2)+((jC(k,2) - jC(j,2))^2)) ;
    if((jC(k,2) - jC(j,2))>=0)
    O = acos((jC(k,1) - jC(j,1))/(L));
    elseif((jC(k,1) - jC(j,1))>=0)
    O=asin((jC(k,2) - jC(j,2))/(L));
    else
    O=pi-(asin((jC(k,2) - jC(j,2))/(L)));    
    end
    T_m=[cos(O) sin(O) 0 0 0 0;-sin(O) cos(O) 0 0 0 0;0 0 1 0 0 0;0 0 0 cos(O) sin(O) 0;0 0 0 (-sin(O)) cos(O) 0;0 0 0 0 0 1] ;
    K_m=((T_m)'*[E(c)*A(c)/L 0 0 -E(c)*A(c)/L 0 0;0 12*(E(c)*I(c)/(L^3)) 6*L*(E(c)*I(c)/(L^3)) 0 -12*(E(c)*I(c)/(L^3)) 6*L*(E(c)*I(c)/(L^3));0 6*L*(E(c)*I(c)/(L^3)) 4*(L^2)*(E(c)*I(c)/(L^3)) 0 -6*L*(E(c)*I(c)/(L^3)) 2*(L^2)*(E(c)*I(c)/(L^3));-E(c)*A(c)/L 0 0 E(c)*A(c)/L 0 0;0 -12*(E(c)*I(c)/(L^3)) -6*L*(E(c)*I(c)/(L^3)) 0 12*(E(c)*I(c)/(L^3)) -6*L*(E(c)*I(c)/(L^3));0 6*L*(E(c)*I(c)/(L^3)) 2*(L^2)*(E(c)*I(c)/(L^3)) 0 -6*L*(E(c)*I(c)/(L^3)) 4*(L^2)*(E(c)*I(c)/(L^3))]*(T_m)) ; %Member Stiffness
    j1=((3*(m_Nj(c,1)))-2); j2=(3*(m_Nj(c,1))-1);j3=(3*(m_Nj(c,1)));  %Index of first joint
    k1=((3*(m_Nj(c,2)))-2); k2=(3*(m_Nj(c,2))-1);k3=(3*(m_Nj(c,2)));  %Index of second joint
    Kj([j1 j2 j3 k1 k2 k3],[j1 j2 j3 k1 k2 k3]) = Kj([j1 j2 j3 k1 k2 k3],[j1 j2 j3 k1 k2 k3])+K_m; %Assembling of Joint Stiffness Matrix 
    A_f(:,c)=(T_m')*[-P_m(c,2)/2 (-P_m(c,1)/2-(P_m(c,3)*L/2)+((P_m(c,4))*3/(2*L))) (-P_m(c,1)*L/8-P_m(c,3)*L*L/12+(P_m(c,4))/4) (-P_m(c,2)/2) (-P_m(c,1)/2-(P_m(c,3)*L/2)-((P_m(c,4))*3/(2*L))) ((P_m(c,1)*L/8)+(P_m(c,3)*L*L/12)+(P_m(c,4))/4)]' ; %Fixed-end actions 
    A_e([j1 j2 j3 k1 k2 k3],1)= A_e([j1 j2 j3 k1 k2 k3],1)-A_f(:,c) ;  
end
A_t = A_j + A_e ; %Total load vector 
%% Rearrangement and Solutions
Kj_re = Kj([FF;BC] ,[FF;BC]) ; %Rearranged Joint Stiffness Matrix
Kf_ff = Kj(FF, FF); % Stiffness matrix corresponding to free DOFS
A_ft = A_t(FF); A_rt = A_t(BC); % Partition of load vector
D_f=(((Kf_ff)^(-1))*(A_ft)); %Displacement at free DOFs 
disp('Displacements at Free DOFs "m"'); disp(D_f);
Kj_rf=Kj_re(length (FF)+1: end , 1:length (FF)); 
R_s = -A_rt + Kj_rf*D_f ;
disp('Support Reactions "N"'); disp (R_s)
D_j = zeros(3*Nj, 1) ;%global displacement vector
D_j(FF) = D_f ;
F_m=zeros (6, m); % Predefining the member end-forces in Global Coordinates system.  
for c=1:m 
    j=m_Nj(c,1) ; k=m_Nj(c,2) ;
    L = sqrt(((jC(k,1) - jC(j,1))^2)+((jC(k,2) - jC(j,2))^2)) ;
    if((jC(k,2) - jC(j,2))>=0)
    O = acos((jC(k,1) - jC(j,1))/(L));
    elseif((jC(k,1) - jC(j,1))>=0)
    O=asin((jC(k,2) - jC(j,2))/(L));
    else
    O=pi-(asin((jC(k,2) - jC(j,2))/(L)));    
    end
    T_m=[cos(O) sin(O) 0 0 0 0;-sin(O) cos(O) 0 0 0 0;0 0 1 0 0 0;0 0 0 cos(O) sin(O) 0;0 0 0 (-sin(O)) cos(O) 0;0 0 0 0 0 1] ;
    K_m=((T_m)'*[E(c)*A(c)/L 0 0 -E(c)*A(c)/L 0 0;0 12*(E(c)*I(c)/(L^3)) 6*L*(E(c)*I(c)/(L^3)) 0 -12*(E(c)*I(c)/(L^3)) 6*L*(E(c)*I(c)/(L^3));0 6*L*(E(c)*I(c)/(L^3)) 4*(L^2)*(E(c)*I(c)/(L^3)) 0 -6*L*(E(c)*I(c)/(L^3)) 2*(L^2)*(E(c)*I(c)/(L^3));-E(c)*A(c)/L 0 0 E(c)*A(c)/L 0 0;0 -12*(E(c)*I(c)/(L^3)) -6*L*(E(c)*I(c)/(L^3)) 0 12*(E(c)*I(c)/(L^3)) -6*L*(E(c)*I(c)/(L^3));0 6*L*(E(c)*I(c)/(L^3)) 2*(L^2)*(E(c)*I(c)/(L^3)) 0 -6*L*(E(c)*I(c)/(L^3)) 4*(L^2)*(E(c)*I(c)/(L^3))]*(T_m)) ;   %Member Stiffness
    j1=((3*(m_Nj(c,1)))-2); j2=(3*(m_Nj(c,1))-1);j3=(3*(m_Nj(c,1)));  %Index of first joint
    k1=((3*(m_Nj(c,2)))-2); k2=(3*(m_Nj(c,2))-1);k3=(3*(m_Nj(c,2)));  %Index of second joint
    D_m(:,1) = D_j([j1 j2 j3 k1 k2 k3],1);
    F_m(:,c)=A_f(:,c)+K_m*D_m;
end
 disp('Member end-forces along the global axis "N"');disp(F_m) 
 L_L=zeros (6, m);% Predefining the member end-forces in Local Coordinates system.  
 for c=1:m 
    j=m_Nj(c,1) ; k=m_Nj(c,2) ;
    L = sqrt(((jC(k,1) - jC(j,1))^2)+((jC(k,2) - jC(j,2))^2)) ;
    if((jC(k,2) - jC(j,2))>=0)
    O = acos((jC(k,1) - jC(j,1))/(L));
    elseif((jC(k,1) - jC(j,1))>=0)
    O=asin((jC(k,2) - jC(j,2))/(L));
    else
    O=pi-(asin((jC(k,2) - jC(j,2))/(L)));    
    end
    T_m=[cos(O) sin(O) 0 0 0 0;-sin(O) cos(O) 0 0 0 0;0 0 1 0 0 0;0 0 0 cos(O) sin(O) 0;0 0 0 (-sin(O)) cos(O) 0;0 0 0 0 0 1] ;
    F_L(:,1) = F_m(:,c); L_L(:,c)=(T_m)*F_L; 
 end
 disp('Member end-forces along the local axis "N"');disp(L_L);
 for c=1:Nj
 j_d(2*c-1) = D_j(3*c-2)+jC(c,1);
 j_d(2*c) = D_j(3*c-1)+jC(c,2);
 end
 Dj = zeros(2*Nj,1); %Joint Coordinates vector  
  for c=1:Nj
 Dj(2*c-1) = jC(c,1);
 Dj(2*c) = jC(c,2);
  end 
for c=1:m
j=m_Nj(c,1) ; k=m_Nj(c,2) ;   
j1=((2*(m_Nj(c,1)))-1); j2=(2*(m_Nj(c,1)));
k1=((2*(m_Nj(c,2)))-1); k2=(2*(m_Nj(c,2)));
L = sqrt(((jC(k,1) - jC(j,1))^2)+((jC(k,2) - jC(j,2))^2));
L1=sqrt(((j_d(k1)-j_d(j1))^2)+((j_d(k2)-j_d(j2))^2));
if c==1
        net__flix__in= jC(j,1);
        net__flix__co__in=jC(j,2);
end
if((jC(k,2) - jC(j,2))>=0)
O = acos((jC(k,1) - jC(j,1))/(L));
elseif((jC(k,1) - jC(j,1))>=0)
O=asin((jC(k,2) - jC(j,2))/(L));
else
O=pi-(asin((jC(k,2) - jC(j,2))/(L)));    
end
if((j_d(k2) -j_d(j2))>=0)
O1 = acos((j_d(k1) -j_d(j1))/(L1));
elseif((j_d(k1) -j_d(j1))>=0) 
O1= asin((j_d(k2) -j_d(j2))/(L1));
else
O1=pi- asin((j_d(k2) -j_d(j2))/(L1));    
end
jx=min(Dj(j1),Dj(k1));kx=max(Dj(j1),Dj(k1));
j__y=min(Dj(j2),Dj(k2));k__y=max(Dj(j2),Dj(k2));
if(kx~=jx)
x = linspace(jx,kx); 
y =(((Dj(j2)-Dj(k2))/(Dj(j1)-Dj(k1)))*(x-(Dj(j1)))+(Dj(j2))) ;
%y=y*scale_factor ;
plot(x,y,":","LineWidth",2.5,"Marker","*","MarkerIndices",[1,(round(length(x)/2)),(length(x))],"MarkerSize",10) ;
else
y = linspace(j__y,k__y); 
x =(((Dj(j1)-Dj(k1))/(Dj(j2)-Dj(k2)))*(y-(Dj(j2)))+(Dj(j1))) ;
%y=y*scale_factor ;
% plot([Dj(j1) Dj(k1)],[Dj(j2) Dj(k2)],":","LineWidth",2.5,"Marker","*","MarkerSize",10);    % ,"MarkerIndices",[1,2,3]
plot(x,y,":","LineWidth",2.5,"Marker","*","MarkerIndices",[1,(round(length(y)/2)),(length(y))],"MarkerSize",10) ;
end    
hold on
x = linspace(0,L/2);
y = ((((L_L(3,c)+L_L(6,c))/L)*(x.^3)/6)-(L_L(3,c)*(x.^2)/2)+(2*L_L(3,c)-L_L(6,c))/6*L*x)/(E(c)*I(c))+((-1*P_m(c,3)/(E(c)*I(c)))*(L*(x.^3)/12-((x.^4)/24)-((L.^3)*x/24)))+((-1*P_m(c,1)/(E(c)*I(c)))*((x.^3)/12-((L.^2)*x/16)))+((P_m(c,4)/(E(c)*I(c)))*((x.^3)/(6*L)-L/24*x)) ;
y=y*scale_factor ;
[muskus,pikachu_h]=max(y);
[muskuss,mickymouse_k]=min(y);
if mickymouse_k == length(y)
    mickymouse_k=1;
end
if pikachu_h == length(y)
    pikachu_h = 1;
end
plot((((x*((L1/2)+(P_m(c,2)*0.25*L/(A(c)*E(c))))/(L/2))*(cos(-O1))+y*(sin(-O1)))+j_d(j1)),((-1*(x*((L1/2)+(P_m(c,2)*0.25*L/(A(c)*E(c))))/(L/2))*(sin(-O1))+y*(cos(-O1)))+j_d(j2)),'b',"LineWidth",2.5,"Marker","*","MarkerIndices",[pikachu_h,mickymouse_k],"MarkerSize",10)  ;
hold on 
x = linspace(L/2,L);
y = ((((L_L(3,c)+L_L(6,c))/L)*(x.^3)/6)-(L_L(3,c)*(x.^2)/2)+(2*L_L(3,c)-L_L(6,c))/6*L*x)/(E(c)*I(c))+((-1*P_m(c,3)/(E(c)*I(c)))*(L*(x.^3)/12-((x.^4)/24)-((L.^3)*x/24)))+((-1*P_m(c,1)/(E(c)*I(c)))*((x.^3)/12-((L.^2)*x/16)-((x-L/2).^3)/6))+((P_m(c,4)/(E(c)*I(c)))*((x.^3)/(6*L)-L/24*x-((x-L/2).^2)/2)) ;
y=y*scale_factor ;
[muskusss,pikachu]=max(y);
[muskussss,mickymouse]=min(y);
if pikachu == 1
    if pikachu_h ~=1
        pikachu=length(y);
    end
end
if mickymouse == 1
    if mickymouse_k ~=1 
        mickymouse=length(y);     
    end
end
plot((((((L1/2)+(P_m(c,2)*0.25*L/(A(c)*E(c))))+(x-L/2)*((L1/2)-(P_m(c,2)*0.25*L/(A(c)*E(c))))/(L/2))*(cos(-O1))+y*(sin(-O1)))+j_d(j1)),((-1*(((L1/2)+(P_m(c,2)*0.25*L/(A(c)*E(c))))+(x-L/2)*((L1/2)-(P_m(c,2)*0.25*L/(A(c)*E(c))))/(L/2))*(sin(-O1))+y*(cos(-O1)))+j_d(j2)),'b',"LineWidth",2.5,"LineStyleMode","auto","Marker","*","MarkerIndices",[pikachu,mickymouse],"MarkerSize",10,'HandleVisibility','off') ;
hold on
% S.F.D
x = linspace(0,L/2);
y =(L_L(2,c))+P_m(c,3)*x;
y=y/scale___factor ;
net_flix_in=(((x*((L/2)+(P_m(c,2)*0*0.25*L/(A(c)*E(c))))/(L/2))*(cos(-O))+y*(sin(-O)))+jC(j,1));
net_flix_co_in=((-1*(x*((L/2)+(P_m(c,2)*0*0.25*L/(A(c)*E(c))))/(L/2))*(sin(-O))+y*(cos(-O)))+jC(j,2));
plot((((x*((L/2)+(P_m(c,2)*0*0.25*L/(A(c)*E(c))))/(L/2))*(cos(-O))+y*(sin(-O)))+jC(j,1)),((-1*(x*((L/2)+(P_m(c,2)*0*0.25*L/(A(c)*E(c))))/(L/2))*(sin(-O))+y*(cos(-O)))+jC(j,2)),'color',[.5 .4 .7],"LineWidth",2.5)  ;
 hold on 
 if c>1
 if m_Nj(c,1)==m_Nj((c-1),2)
   plot([net__flix__in(length(net__flix__in)) net_flix_in((1))],[net__flix__co__in(length(net__flix__co__in)) net_flix_co_in((1))],'color',[.5 .4 .7],"LineWidth",2.5,'HandleVisibility','off');
 else
    plot([jC(j,1) net_flix_in((1))],[jC(j,2) net_flix_co_in((1))],'color',[.5 .4 .7],"LineWidth",2.5,'HandleVisibility','off'); 
    plot([jC(m_Nj((c-1),2),1) net__flix__in(length(net__flix__in))],[jC(m_Nj((c-1),2),2) net__flix__co__in(length(net__flix__co__in))],'color',[.5 .4 .7],"LineWidth",2.5,'HandleVisibility','off');
 end
 else
   plot([net__flix__in(length(net__flix__in)) net_flix_in((1))],[net__flix__co__in(length(net__flix__co__in)) net_flix_co_in((1))],'color',[.5 .4 .7],"LineWidth",2.5,'HandleVisibility','off');  
 end
hold on
x = linspace(L/2,L);
y= (L_L(2,c))+P_m(c,3)*x+P_m(c,1);
y=y/scale___factor ;
net__flix__in= (((((L/2)+(P_m(c,2)*0*0.25*L/(A(c)*E(c))))+(x-L/2)*((L/2)-(P_m(c,2)*0*0.25*L/(A(c)*E(c))))/(L/2))*(cos(-O))+y*(sin(-O)))+jC(j,1)) ;
net__flix__co__in = ((-1*(((L/2)+(P_m(c,2)*0*0.25*L/(A(c)*E(c))))+(x-L/2)*((L/2)-(P_m(c,2)*0*0.25*L/(A(c)*E(c))))/(L/2))*(sin(-O))+y*(cos(-O)))+jC(j,2))  ;
plot((((((L/2)+(P_m(c,2)*0*0.25*L/(A(c)*E(c))))+(x-L/2)*((L/2)-(P_m(c,2)*0*0.25*L/(A(c)*E(c))))/(L/2))*(cos(-O))+y*(sin(-O)))+jC(j,1)),((-1*(((L/2)+(P_m(c,2)*0*0.25*L/(A(c)*E(c))))+(x-L/2)*((L/2)-(P_m(c,2)*0*0.25*L/(A(c)*E(c))))/(L/2))*(sin(-O))+y*(cos(-O)))+jC(j,2)),'color',[.5 .4 .7],"LineWidth",2.5,"LineStyleMode","auto",'HandleVisibility','off') ;
hold on 
     plot([net__flix__in(1) net_flix_in(length(net_flix_in))],[net__flix__co__in(1) net_flix_co_in(length(net_flix_co_in))],'color',[.5 .4 .7],"LineWidth",2.5,'HandleVisibility','off');
hold on 
if c==m
    net_flix_in = jC(k,1);
    net_flix_co_in =jC(k,2);
    plot([net__flix__in(length(net__flix__in)) net_flix_in((1))],[net__flix__co__in(length(net__flix__co__in)) net_flix_co_in((1))],'color',[.5 .4 .7],"LineWidth",2.5,'HandleVisibility','off');
hold on
end
% B.M.D
if c==1
        net__flix__in_bmd= jC(j,1);
        net__flix__co__in_bmd=jC(j,2);
end
x = linspace(0,L/2);
y =(-1*(L_L(3,c))+((L_L(6,c))+(L_L(3,c)))/L*x)-P_m(c,1)/2*x-(P_m(c,3)*(L*x/2-x.*x/2))+(P_m(c,4)/L*x);
y=y/scale_____factor ;
net_flix_in_bmd=(((x*((L/2)+(P_m(c,2)*0*0.25*L/(A(c)*E(c))))/(L/2))*(cos(-O))+y*(sin(-O)))+jC(j,1));
net_flix_co_in_bmd=((-1*(x*((L/2)+(P_m(c,2)*0*0.25*L/(A(c)*E(c))))/(L/2))*(sin(-O))+y*(cos(-O)))+jC(j,2));
plot((((x*((L/2)+(P_m(c,2)*0*0.25*L/(A(c)*E(c))))/(L/2))*(cos(-O))+y*(sin(-O)))+jC(j,1)),((-1*(x*((L/2)+(P_m(c,2)*0*0.25*L/(A(c)*E(c))))/(L/2))*(sin(-O))+y*(cos(-O)))+jC(j,2)),'color','g',"LineWidth",2.5)  ;
hold on
 if c>1
 if m_Nj(c,1)==m_Nj((c-1),2)
plot([net__flix__in_bmd(length(net__flix__in_bmd)) net_flix_in_bmd((1))],[net__flix__co__in_bmd(length(net__flix__co__in_bmd)) net_flix_co_in_bmd((1))],'color','g',"LineWidth",2.5,'HandleVisibility','off');
hold on
 else
    plot([jC(j,1) net_flix_in_bmd((1))],[jC(j,2) net_flix_co_in_bmd((1))],'color','g',"LineWidth",2.5,'HandleVisibility','off'); 
    hold on 
     plot([jC(m_Nj((c-1),2),1) net__flix__in_bmd(length(net__flix__in_bmd))],[jC(m_Nj((c-1),2),2) net__flix__co__in_bmd(length(net__flix__co__in_bmd))],'color','g',"LineWidth",2.5,'HandleVisibility','off');
    hold on
 end
 else
     plot([net__flix__in_bmd(length(net__flix__in_bmd)) net_flix_in_bmd((1))],[net__flix__co__in_bmd(length(net__flix__co__in_bmd)) net_flix_co_in_bmd((1))],'color','g',"LineWidth",2.5,'HandleVisibility','off');
hold on
 end 
x = linspace(L/2,L);
y = (-1*(L_L(3,c))+((L_L(6,c))+(L_L(3,c)))/L*x)-(P_m(c,1)/2*(L-x))-(P_m(c,3)*(L*x/2-x.*x/2))+(P_m(c,4)*(x/L-1));
y=y/scale_____factor ;
net__flix__in_bmd= (((((L/2)+(P_m(c,2)*0*0.25*L/(A(c)*E(c))))+(x-L/2)*((L/2)-(P_m(c,2)*0*0.25*L/(A(c)*E(c))))/(L/2))*(cos(-O))+y*(sin(-O)))+jC(j,1)) ;  
net__flix__co__in_bmd =  ((-1*(((L/2)+(P_m(c,2)*0*0.25*L/(A(c)*E(c))))+(x-L/2)*((L/2)-(P_m(c,2)*0*0.25*L/(A(c)*E(c))))/(L/2))*(sin(-O))+y*(cos(-O)))+jC(j,2))   ;
plot((((((L/2)+(P_m(c,2)*0*0.25*L/(A(c)*E(c))))+(x-L/2)*((L/2)-(P_m(c,2)*0*0.25*L/(A(c)*E(c))))/(L/2))*(cos(-O))+y*(sin(-O)))+jC(j,1)),((-1*(((L/2)+(P_m(c,2)*0*0.25*L/(A(c)*E(c))))+(x-L/2)*((L/2)-(P_m(c,2)*0*0.25*L/(A(c)*E(c))))/(L/2))*(sin(-O))+y*(cos(-O)))+jC(j,2)),'color','g',"LineWidth",2.5,"LineStyleMode","auto",'HandleVisibility','off') ;
hold on 
plot([net__flix__in_bmd(1) net_flix_in_bmd(length(net_flix_in_bmd))],[net__flix__co__in_bmd(1) net_flix_co_in_bmd(length(net_flix_co_in_bmd))],'color','g',"LineWidth",2.5,'HandleVisibility','off');
hold on 
if c==m
    net_flix_in_bmd = jC(k,1);
    net_flix_co_in_bmd =jC(k,2);
    plot([net__flix__in_bmd(length(net__flix__in_bmd)) net_flix_in_bmd((1))],[net__flix__co__in_bmd(length(net__flix__co__in_bmd)) net_flix_co_in_bmd((1))],'color','g',"LineWidth",2.5,'HandleVisibility','off');
hold on
end
end
hold off
% axis equal 
xlabel("length")
ylabel("Height")
str1="deflected shape in meter, Scale factor = ";
str__1="S.F.D in 'N', Scale factor = " ;
str___1="B.M.D in 'N-m', Scale factor = " ;
str2=num2str(scale_factor);
str__2=num2str((1/scale___factor)); 
str___2=num2str((1/scale_____factor)); 
str = append(str1,str2);
str__str = append(str__1,str__2);
str___str = append(str___1,str___2);
title("Plot of Original shape with it's deflected shape")
legend("Original shape",str,str__str,str___str)   
n=100;
x_y = zeros(2,n+1);
for c=1:m 
j=m_Nj(c,1) ; k=m_Nj(c,2) ;   
j1=((2*(m_Nj(c,1)))-1); j2=(2*(m_Nj(c,1)));
k1=((2*(m_Nj(c,2)))-1); k2=(2*(m_Nj(c,2)));
L = sqrt(((jC(k,1) - jC(j,1))^2)+((jC(k,2) - jC(j,2))^2));
L1=sqrt(((j_d(k1)-j_d(j1))^2)+((j_d(k2)-j_d(j2))^2));
if((jC(k,2) - jC(j,2))>=0)
O = acos((jC(k,1) - jC(j,1))/(L));
elseif((jC(k,1) - jC(j,1))>=0)
O=asin((jC(k,2) - jC(j,2))/(L));
else
O=pi-(asin((jC(k,2) - jC(j,2))/(L)));    
end
if((j_d(k2) -j_d(j2))>=0)
O1 = acos((j_d(k1) -j_d(j1))/(L1));
elseif((j_d(k1) -j_d(j1))>=0)
O1= asin((j_d(k2) -j_d(j2))/(L1));
else
O1=pi- asin((j_d(k2) -j_d(j2))/(L1));    
end
ds=L/n;
f=(0:ds:L)' ;
g = zeros(n+1,1);
for i=1:n+1
    if f(i)>=0 && f(i)<(L/2)
     g(i) = ((((L_L(3,c)+L_L(6,c))/L)*(f(i).^3)/6)-(L_L(3,c)*(f(i).^2)/2)+(2*L_L(3,c)-L_L(6,c))/6*L*f(i))/(E(c)*I(c))+((-1*P_m(c,3)/(E(c)*I(c)))*(L*(f(i).^3)/12-((f(i).^4)/24)-((L.^3)*f(i)/24)))+((-1*P_m(c,1)/(E(c)*I(c)))*((f(i).^3)/12-((L.^2)*f(i)/16)))+((P_m(c,4)/(E(c)*I(c)))*((f(i).^3)/(6*L)-L/24*f(i)));
     xx=(((f(i)*((L1/2)+(P_m(c,2)*0.25*L/(A(c)*E(c))))/(L/2))*(cos(-O1))+g(i)*(sin(-O1)))+j_d(j1));
     yy=((-1*(f(i)*((L1/2)+(P_m(c,2)*0.25*L/(A(c)*E(c))))/(L/2))*(sin(-O1))+g(i)*(cos(-O1)))+j_d(j2));
     x_y(1,i)=xx;x_y(2,i)=yy;
    else
     g(i)= ((((L_L(3,c)+L_L(6,c))/L)*(f(i).^3)/6)-(L_L(3,c)*(f(i).^2)/2)+(2*L_L(3,c)-L_L(6,c))/6*L*f(i))/(E(c)*I(c))+((-1*P_m(c,3)/(E(c)*I(c)))*(L*(f(i).^3)/12-((f(i).^4)/24)-((L.^3)*f(i)/24)))+((-1*P_m(c,1)/(E(c)*I(c)))*((f(i).^3)/12-((L.^2)*f(i)/16)-((f(i)-L/2).^3)/6))+((P_m(c,4)/(E(c)*I(c)))*((f(i).^3)/(6*L)-L/24*f(i)-((f(i)-L/2).^2)/2));   
     xx=(((((L1/2)+(P_m(c,2)*0.25*L/(A(c)*E(c))))+(f(i)-L/2)*((L1/2)-(P_m(c,2)*0.25*L/(A(c)*E(c))))/(L/2))*(cos(-O1))+g(i)*(sin(-O1)))+j_d(j1));
     yy=((-1*(((L1/2)+(P_m(c,2)*0.25*L/(A(c)*E(c))))+(f(i)-L/2)*((L1/2)-(P_m(c,2)*0.25*L/(A(c)*E(c))))/(L/2))*(sin(-O1))+g(i)*(cos(-O1)))+j_d(j2))  ;
     x_y(1,i)=xx;x_y(2,i)=yy;       
    end
end
disp('Structure Displaced Joint Coordinates=======>>>(x(m),y(m))');disp(c);disp(x_y);
end
n=100;
x_y = zeros(2,n+1);
for c=1:m
j=m_Nj(c,1) ; k=m_Nj(c,2) ;   
j1=((2*(m_Nj(c,1)))-1); j2=(2*(m_Nj(c,1)));
k1=((2*(m_Nj(c,2)))-1); k2=(2*(m_Nj(c,2)));
L = sqrt(((jC(k,1) - jC(j,1))^2)+((jC(k,2) - jC(j,2))^2));
L1=sqrt(((j_d(k1)-j_d(j1))^2)+((j_d(k2)-j_d(j2))^2));
if((jC(k,2) - jC(j,2))>=0)
O = acos((jC(k,1) - jC(j,1))/(L));
elseif((jC(k,1) - jC(j,1))>=0)
O=asin((jC(k,2) - jC(j,2))/(L));
else
O=pi-(asin((jC(k,2) - jC(j,2))/(L)));    
end
if((j_d(k2) -j_d(j2))>=0)
O1 = acos((j_d(k1) -j_d(j1))/(L1));
elseif((j_d(k1) -j_d(j1))>=0)
O1= asin((j_d(k2) -j_d(j2))/(L1));
else
O1=pi- asin((j_d(k2) -j_d(j2))/(L1));    
end
jx=min(Dj(j1),Dj(k1));kx=max(Dj(j1),Dj(k1));
jy=min(Dj(j2),Dj(k2));ky=max(Dj(j2),Dj(k2));
ds=L/n;
f=(0:ds:L)' ;
g = zeros(n+1,1);
for i=1:n+1
    if(kx~=jx)
    g(i) =(((Dj(j2)-Dj(k2))/(Dj(j1)-Dj(k1)))*(f(i)-(Dj(j1)))+(Dj(j2)));
     x_y(1,i)= f(i); x_y(2,i)= g(i);
    else
     x_y(1,i)=kx;x_y(2,i)=jy+(i-1)*ds;       
    end
end
disp('Structure Joint Coordinates=======>>>(x(m),y(m))');disp(c);disp(x_y);
end
%%