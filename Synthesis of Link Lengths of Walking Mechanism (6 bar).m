clc;clear;close;

%Initialisation

p1 = -4.84-6.37i;  % at theta2 = 0 degree
p2 = -7.75-6.6i;   % at theta2 = 80 degree
p3 = -12.39-7.09i; % at theta2 = 180 degree
p4 = -11.02-4.68i; % at theta2 = 225 degree
p5 = -8.08-4.36i;  % at theta2 = 280 degree

theta2_01 = 80*pi/180;
theta2_02 = 180*pi/180;
theta2_03 = 225*pi/180;
theta2_04 = 280*pi/180;

l1 = 5*exp(15*pi/180*1i); % Assuming l1
%% Solving loop-closure 1
syms l2 l3 theta3_01 theta3_02 theta3_03   % Declaring the unknowns

eq1 = l1*conj(l2)*(exp(-theta2_01*1i)-1)+l1*conj(l3)*(exp(-theta3_01*1i)-1)+l2*conj(l1)*(exp(theta2_01*1i)-1)+l3*conj(l1)*(exp(theta3_01*1i)-1)+l2*conj(l3)*(exp((theta2_01-theta3_01)*1i)-1)+ l3*conj(l2)*(exp((theta3_01-theta2_01)*1i)-1); 
eq2 = l1*conj(l2)*(exp(-theta2_02*1i)-1)+l1*conj(l3)*(exp(-theta3_02*1i)-1)+l2*conj(l1)*(exp(theta2_02*1i)-1)+l3*conj(l1)*(exp(theta3_02*1i)-1)+l2*conj(l3)*(exp((theta2_02-theta3_02)*1i)-1)+ l3*conj(l2)*(exp((theta3_02-theta2_02)*1i)-1); 
eq3 = l1*conj(l2)*(exp(-theta2_03*1i)-1)+l1*conj(l3)*(exp(-theta3_03*1i)-1)+l2*conj(l1)*(exp(theta2_03*1i)-1)+l3*conj(l1)*(exp(theta3_03*1i)-1)+l2*conj(l3)*(exp((theta2_03-theta3_03)*1i)-1)+ l3*conj(l2)*(exp((theta3_03-theta2_03)*1i)-1); 

sol = vpasolve([eq1 eq2 eq3],[l2 l3 theta3_01 theta3_02 theta3_03],[1.8 -5+0.2i 0.21 0.09 -0.35]);

l4 = (l1+sol.l2+sol.l3);
l2 = (sol.l2*exp(0*1i));
l3 = (sol.l3);

%[abs(l2) abs(l3) abs(l4)]
%% Solving loop closure 2
c1 = p1-l2-l3;
c2 = p2-l2*exp(theta2_01*1i)-l3*exp(sol.theta3_01*1i);
c3 = p3-l2*exp(theta2_02*1i)-l3*exp(sol.theta3_02*1i);

syms a3 w theta6_01 theta6_02 alpha beta % Declaring the unknowns

eq4 = a3*cos(alpha)+w*cos(beta)==real(c1);
eq5 = a3*sin(alpha)+w*sin(beta)==imag(c1);
eq6 = a3*cos(sol.theta3_01+alpha)+w*cos(beta+theta6_01)==real(c2);
eq7 = a3*sin(sol.theta3_01+alpha)+w*sin(beta+theta6_01)==imag(c2);
eq8 = a3*cos(sol.theta3_02+alpha)+w*cos(beta+theta6_02)==real(c3);
eq9 = a3*sin(sol.theta3_02+alpha)+w*sin(beta+theta6_02)==imag(c3);


sol1 = vpasolve([eq4 eq5 eq6 eq7 eq8 eq9],[a3 w theta6_01 theta6_02 alpha beta],[3.5 7.5 -0.2 -0.5 2.93 4.939]);
%% Finding unknowns for loop-closure 3
syms theta3_04 theta6_03 theta6_04 % Declaring the unknowns

theta3_04 = vpasolve(l1*conj(l2)*(exp(-theta2_04*1i)-1)+l1*conj(l3)*(exp(-theta3_04*1i)-1)+l2*conj(l1)*(exp(theta2_04*1i)-1)+l3*conj(l1)*(exp(theta3_04*1i)-1)+l2*conj(l3)*(exp((theta2_04-theta3_04)*1i)-1)+ l3*conj(l2)*(exp((theta3_04-theta2_04)*1i)-1),theta3_04,-0.5);

theta4_01 = acos(real((l1+l2*exp(theta2_01*1i)+l3*exp(sol.theta3_01*1i))/(l1+l2+l3)));
theta4_02 = acos(real((l1+l2*exp(theta2_02*1i)+l3*exp(sol.theta3_02*1i))/(l1+l2+l3)));
theta4_03 = acos(real(vpa((l1+l2*exp(theta2_03*1i)+l3*exp(sol.theta3_03*1i))/(l1+l2+l3))));
theta4_04 = acos(real((l1+l2*exp(theta2_04*1i)+l3*exp(theta3_04*1i))/(l1+l2+l3)));

c4 = vpa(p4-l2*exp(theta2_03*1i)-l3*exp(sol.theta3_03*1i));
c5 = p5-l2*exp(theta2_04*1i)-l3*exp(theta3_04*1i);


theta6_03 = vpasolve(sol1.a3*cos(sol.theta3_03+sol1.alpha)+sol1.w*cos(sol1.beta+theta6_03)==real(c4),theta6_03,-0.53);
theta6_04 = vpasolve(sol1.a3*cos(theta3_04+sol1.alpha)+sol1.w*cos(sol1.beta+theta6_04)==real(c5),theta6_04,-0.43);
%% Solving loop-closure 3
syms l5 l6 theta5_01 theta5_02 theta5_03 theta5_04 % Declaring the unknowns

eq10 = l6*(exp(sol1.theta6_01)-1)+l5*(exp(theta5_01)-1) == l4*(exp(theta4_01)-1)+sol1.a3*exp(sol1.alpha*1i)*(exp(sol.theta3_01)-1);
eq11 = l6*(exp(sol1.theta6_02)-1)+l5*(exp(theta5_02)-1) == l4*(exp(theta4_02)-1)+sol1.a3*exp(sol1.alpha*1i)*(exp(sol.theta3_01)-1);
eq12 = l6*(exp(theta6_03)-1)+l5*(exp(theta5_03)-1) == l4*(exp(theta4_03)-1)+sol1.a3*exp(sol1.alpha*1i)*(exp(sol.theta3_01)-1);
eq13 = l6*(exp(theta6_04)-1)+l5*(exp(theta5_04)-1) == l4*(exp(theta4_04)-1)+sol1.a3*exp(sol1.alpha*1i)*(exp(theta3_04)-1);

sol2 = vpasolve([eq10 eq11 eq12 eq13],[l5 l6 theta5_01 theta5_02 theta5_03 theta5_04],[-2.4+3.6i -2.6-4.2i 0.2 0.76 0.12 -0.37]);

a1 = (l4+sol1.a3*exp(sol1.alpha*1i)-sol2.l6-sol2.l5);

%% Displaying the link lengths

disp(['l1 = ', num2str(l1),' ,  |l1| = ',num2str(abs(l1))]);
disp(['l2 = ', num2str(double(l2)),' ,  |l2| = ',num2str(abs(double(l2)))]);
disp(['l3 = ', num2str(double(l3)),' ,  |l3| = ',num2str(abs(double(l3)))]);
disp(['l4 = ', num2str(double(l4)),' ,  |l4| = ',num2str(abs(double(l4)))]);
disp(['l5 = ', num2str(double(sol2.l5)),' ,  |l5| = ',num2str(abs(double(sol2.l5)))]);
disp(['l6 = ', num2str(double(sol2.l6)),' ,  |l6| = ',num2str(abs(double(sol2.l6)))]);
disp(['a3l3 = ', num2str(double(sol1.a3*exp(sol1.alpha*1i))),' ,  |a3l3| = ',num2str(abs(double(sol1.a3*exp(sol1.alpha*1i))))]);
disp(['a6l6 = ', num2str(double(sol1.w*exp(sol1.beta*1i))),' ,  |a6l6| = ',num2str(abs(double(sol1.w*exp(sol1.beta*1i))))]);










