close all
clear all

load NACA0012_flowfieldv2.dat

imax = 513;
jmax = 257;
TE_start = 65; % Trailing Edge Lower Point
TE_end = 449; % Trailing Edge Upper Point
LE = 257; % Leading Edge Point
gamma = 1.4;
mach = 0.1;
p0 = 1.; % Freestream static pressure
alpha = 8; % Angle of Attack
Re = 3e6; % Reynolds number



k = 1;
for j=1:jmax
    for i=1:imax
        x(i,j) = NACA0012_flowfieldv2(k,1);
        y(i,j) = NACA0012_flowfieldv2(k,2);
        k = k +1;
    end
end



% Plot the Airfoil Surface
figure(1)
plot(x(TE_start:TE_end,1),...
    y(TE_start:TE_end,1),'k*-')
axis equal



%Plot the Computational Domain.
figure(2)
axis equal
hold on
for j=1:jmax
    plot(x(1:imax,j),y(1:imax,j),'k-')  
end
for i=1:imax
    plot(x(i,1:jmax),y(i,1:jmax),'k-')
end


%Plot the Computational Domain.
figure(3)
axis equal
hold on
for j=1:jmax
    plot(x(1:imax,j),y(1:imax,j),'k-')  
end
for i=1:imax
    plot(x(i,1:jmax),y(i,1:jmax),'k-')
end
axis([-.1 1.1 -.2 .2])


Nvariables = 7;
% Transfer all other flow properties and evaluate Pressure
k = 1;
for j=1:jmax
    for i=1:imax
        for n=1:Nvariables
            w(i,j,n) = NACA0012_flowfieldv2(k,2+n);
        end
        k = k +1;
    end
end




% ***********************************%
% Plot Surface Pressure
% ***********************************%
for j=1:jmax
    for i=1:imax
        pressure(i,j) = (gamma -1)*(w(i,j,4) ...
            -.5*(w(i,j,2)^2  +w(i,j,3)^2)/w(i,j,1));
        velocity(i,j) = sqrt((w(i,j,2)/w(i,j,1))^2 ...
                            +(w(i,j,3)/w(i,j,1))^2);
    end 
end


% Evaluate Shear Stress.
for i=TE_start:TE_end
    dxi = x(i,1) -x(i-1,1); % dx/dxi
    dyi = y(i,1) -y(i-1,1); % dy/dxi
    dxj = x(i,2) -x(i,1); % dx/deta
    dyj = y(i,2) -y(i,1); % dy/deta
    dsj = 1/(dxi*dyj -dyi*dxj);

    dui = w(i,1,2)/w(i,1,1) -w(i-1,1,2)/w(i-1,1,1); % du/dxi
    dvi = w(i,1,3)/w(i,1,1) -w(i-1,1,3)/w(i-1,1,1); % dv/dxi 
    duj = w(i,2,2)/w(i,2,1) -w(i,1,2)/w(i,1,1); % du/deta
    dvj = w(i,2,3)/w(i,2,1) -w(i,1,3)/w(i,1,1); % dv/deta

    dux = (dui*dyj -duj*dyi)*dsj; % du/dx
end








