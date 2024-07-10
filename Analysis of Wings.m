clc;clear;close;
l1 =79.204 ; l2 =54.00 ; l3 =57.17 ; l4 = 24.696;
a =1.5; alpha =pi;
theta1 =4.2214 ; N =100;

%% Solving nonlinear equations, storing the solutions, and updating initial guess
theta_2_and_3 = [pi,pi];
theta4 = linspace(0,4*pi,N);
theta2 = 0*theta4;
theta3 = 0*theta4;
for i=1:N
    xsol = fsolve(@(x)loop_closure_kdm_wing(x,l1,l2,l3,l4,theta1,theta4(i)),theta_2_and_3);
    theta2(i) = xsol(1);
    theta3(i) = xsol(2);
    theta_2_and_3 = [ theta2(i), theta3(i)];
   
end
coupler_curve_x_coord = l2*cos(theta2);
coupler_curve_y_coord = l2*sin(theta2);
P_x_coord_1 = a*l2*cos(theta2-alpha);
P_y_coord_1 = a*l2*sin(theta2-alpha);

%% Plotting the motion of the mechanism and draw coupler curve  
for i=1:N
    x_coord = [0,l2*cos(theta2(i)),l2*cos(theta2(i))+l3*cos(theta3(i)),l2*cos(theta2(i))+l3*cos(theta3(i))+l4*cos(theta4(i))];
    y_coord = [0,l2*sin(theta2(i)),l2*sin(theta2(i))+l3*sin(theta3(i)),l2*sin(theta2(i))+l3*sin(theta3(i))+l4*sin(theta4(i))];
    P_x_coord = [0,a*l2*cos(theta2(i)-alpha)];
    P_y_coord = [0,a*l2*sin(theta2(i)-alpha)];
    figure(3)
    plot(x_coord,y_coord, 'k-',LineWidth=4)
    hold on;
    plot(P_x_coord,P_y_coord, 'b--',LineWidth=3)
    %plot(coupler_curve_x_coord(1:i),coupler_curve_y_coord(1:i), 'r-')
    plot(P_x_coord_1(1:i),P_y_coord_1(1:i),'r-', LineWidth=6)
    plot([x_coord,P_x_coord],[y_coord,P_y_coord], 'ro')
    hold off;
    axis([-l1-l4 2*l1+l4 -l1-l4 l1+l4])
    grid on
    title("kdm wing",'FontSize',30,'FontName','Palatino Linotype')
    xlabel('$x$','interpreter','latex','FontSize',30,'FontName','Palatino Linotype')
    ylabel('$y$','interpreter','latex','FontSize',30,'FontName','Palatino Linotype')
    pause(0.1)
end

%% Loop-closure equations defined as a function for fsolve() to call
function F=loop_closure_kdm_wing(x,l1,l2,l3,l4,theta1,theta4)
F(1) = l2*cos(x(1))+l3*cos(x(2))+l4*cos(theta4)-l1*cos(theta1);
F(2) = l2*sin(x(1))+l3*sin(x(2))+l4*sin(theta4)-l1*sin(theta1);
end