%% 90 x 80 Grid (Medium Grid) 

%% clear
% clear
% clc

%%
%X-DIRECTION
 %Constant grid spacing in x-direction above airfoil
  Nx_2=50; 
  x_coord2=linspace(20, 21, Nx_2);
  dx_2=1/(Nx_2-1);

   Nx_1=20; 
   dx_1=.35;
   x_1=20;
   x_coord1=zeros(1,Nx_1);
   for i=1:Nx_1
       dx_1=exp(0.01*i)*dx_1;
       x_coord1(Nx_1-(i-1))=x_1-dx_1;
       x_1=x_coord1(Nx_1-(i-1));
   end


 Nx_3=20; 
   dx_3=1.45;
   x_3=21;
   x_coord3=zeros(1,Nx_3);
     x_coord3=zeros(1,Nx_3);
   for i=1:Nx_3
       dx_3=exp(0.000009*i)*dx_3;
       x_coord3(i)=x_3+dx_3;
       x_3=x_coord3(i);
   end
   
   Nx=Nx_1+Nx_2+Nx_3;
   x=[x_coord1, x_coord2, x_coord3];
   
   %Y-DIRECTION
   %Exponential stretching of the grid in y-direction
   Ny=80; 
   
   tc = 0.08;
   dy=tc/2;
   y=zeros(1, 20); 
   y(2)=dy;
   for i=3:Ny
       dy=exp(0.0014862*i)*dy;
       y(i)=y(i-1)+dy;
   end
   
dy=tc/2;

%% Specifying Parameters


x_length = length(x);
y_length = length(y);
phi = zeros(x_length,y_length);
A = zeros(x_length,y_length);
mu = zeros(x_length,y_length);
a = zeros(x_length,y_length);
b = zeros(x_length,y_length);
c = zeros(x_length,y_length);
d = zeros(x_length,y_length);
e = zeros(x_length,y_length);
g = zeros(x_length,y_length);
LE = 21;
TE = 70;
M = 0.80;
gamma = 1.4;
R = 287.058;
T = 293;
U = M * sqrt(gamma*R*T);

%% Initialising Boundary Conditions

phi(1,:) = 0;
phi(x_length,:) = 0;
phi(:,y_length) = 0;
phi(1:LE-1,1) = phi(1:LE-1,2);
phi(TE+1:x_length,1) = phi(TE+1:x_length,2);
le = LE;
te = TE;
omega = 0.5;
for i=21:70
    y1(i) = 0.08*(-2*x(i)^2 + 82*x(i) - 840);% + omega*DELS80(i-20);
    
end
y1(71) = 0;
for l=21:70
    
%         if l == 21
%             dydx = (y1(l+1) - y1(l))/(x(l+1) - x(l));
%         elseif l==70
%             dydx = (y1(l) - y1(l-1))/(x(l) - x(l-1));
%         else
            dydx = (y1(l+1) - y1(l-1))/(x(l+1) - x(l-1));
%         end
    phi(l,1) = -dy*U*dydx + phi(l,2);
  
end

for i=3:x_length-1
    for j=2:y_length-1
        
        A(i,j) = (1 - M^2) - ((gamma + 1)*(M^2)*(phi(i+1,j) - phi(i-1,j)))/(U*(x(i+1) - x(i-1)));
        
        if A(i,j)>0
            mu(i,j) = 0;
        end
        if A(i,j)<0
            mu(i,j) = 1;
        end
    end
end

%% First loop for phi

for i = 3:x_length-1
    for j = 2:y_length-1
        
        a(i,j) = (-2*(1 - mu(i,j))*A(i,j))/((x(i+1) - x(i))*(x(i+1) - x(i-1))) + (-2*(1 - mu(i,j))*A(i,j))/((x(i) - x(i-1))*(x(i+1) - x(i-1))) + (2*mu(i-1,j)*A(i-1,j))/((x(i) - x(i-1))*(x(i) - x(i-2))) + (-2)/((y(j+1) - y(j))*(y(j+1) - y(j-1))) + (-2)/((y(j) - y(j-1))*(y(j+1) - y(j-1)));
        b(i,j) = 2/((y(j+1) - y(j))*(y(j+1) - y(j-1)));
        c(i,j) = 2/((y(j) - y(j-1))*(y(j+1) - y(j-1)));
        d(i,j) = (2*(1 - mu(i,j))*A(i,j))/((x(i) - x(i-1))*(x(i+1) - x(i-1))) + (-2*mu(i-1,j)*A(i-1,j))/((x(i) - x(i-1))*(x(i) - x(i-2))) + (-2*mu(i-1,j)*A(i-1,j))/((x(i-1) - x(i-2))*(x(i) - x(i-2)));
        e(i,j) = (2*(1 - mu(i,j))*A(i,j))/((x(i+1) - x(i))*(x(i+1) - x(i-1)));
        g(i,j) = (2*mu(i-1,j)*A(i-1,j))/((x(i-1) - x(i-2))*(x(i) - x(i-2)));
        
        phi(i,j) = (-c(i,j)*phi(i,j-1) - g(i,j)*phi(i-2,j) - d(i,j)*phi(i-1,j) - e(i,j)*phi(i+1,j) - b(i,j)*phi(i,j+1))/a(i,j);
        
    end
end
phiITER = phi;

%% Gauss-Siedel Method

err = 1;
iterSiedel = 1;
tic
while err > 1*10^-4
    
    phi(1,:) = 0;
    phi(2,:) = 0;
    phi(x_length,:) = 0;
    phi(:,y_length) = 0;
    phi(1:LE-1,1) = phi(1:LE-1,2);
    phi(TE+1:x_length,1) = phi(TE+1:x_length,2);
    le = LE;
    te = TE;

    for l=21:70
%         if l == 21
%             dydx = (y1(l+1) - y1(l))/(x(l+1) - x(l));
%         elseif l==70
%             dydx = (y1(l) - y1(l-1))/(x(l) - x(l-1));
%         else
            dydx = (y1(l+1) - y1(l-1))/(x(l+1) - x(l-1));
%         end
            phi(l,1) = -dy*U*dydx + phi(l,2);

    end

    for i=3:x_length-1
        for j=2:y_length-1

            A(i,j) = (1 - M^2) - ((gamma + 1)*(M^2)*(phi(i+1,j) - phi(i-1,j)))/(U*(x(i+1) - x(i-1)));

            if A(i,j)>0
                mu(i,j) = 0;
            end
            if A(i,j)<0
                mu(i,j) = 1;
            end
        end
    end
    
    for i = 3:x_length-1
        for j = 2:y_length-1

            a(i,j) = (-2*(1 - mu(i,j))*A(i,j))/((x(i+1) - x(i))*(x(i+1) - x(i-1))) + (-2*(1 - mu(i,j))*A(i,j))/((x(i) - x(i-1))*(x(i+1) - x(i-1))) + (2*mu(i-1,j)*A(i-1,j))/((x(i) - x(i-1))*(x(i) - x(i-2))) + (-2)/((y(j+1) - y(j))*(y(j+1) - y(j-1))) + (-2)/((y(j) - y(j-1))*(y(j+1) - y(j-1)));
            b(i,j) = 2/((y(j+1) - y(j))*(y(j+1) - y(j-1)));
            c(i,j) = 2/((y(j) - y(j-1))*(y(j+1) - y(j-1)));
            d(i,j) = (2*(1 - mu(i,j))*A(i,j))/((x(i) - x(i-1))*(x(i+1) - x(i-1))) + (-2*mu(i-1,j)*A(i-1,j))/((x(i) - x(i-1))*(x(i) - x(i-2))) + (-2*mu(i-1,j)*A(i-1,j))/((x(i-1) - x(i-2))*(x(i) - x(i-2)));
            e(i,j) = (2*(1 - mu(i,j))*A(i,j))/((x(i+1) - x(i))*(x(i+1) - x(i-1)));
            g(i,j) = (2*mu(i-1,j)*A(i-1,j))/((x(i-1) - x(i-2))*(x(i) - x(i-2)));

            phi(i,j) = (-c(i,j)*phi(i,j-1) - g(i,j)*phi(i-2,j) - d(i,j)*phi(i-1,j) - e(i,j)*phi(i+1,j) - b(i,j)*phi(i,j+1))/a(i,j);

        end
    end
    
    err = max(max(abs(phiITER - phi)))
    residualSiedel(iterSiedel) = max(max(abs(phiITER - phi)));
    phiITER = phi;
    iterSiedel = iterSiedel + 1;
end

% contour(x,y,transpose(phi)); hold on; contour(x,-y,transpose(phi)); hold off;
% xlabel('X - axis \rightarrow');
% ylabel('Y - axis \rightarrow');
% title('Contour Plot of \phi (M_\infty = 0.88)');

z = 1;
for i = 21:70
       j = 1;
          
          cp80(z) = -2*(phi(i,j) - phi(i-1,j))/((x(i) - x(i-1))*U);
          u(z) = (phi(i,j) - phi(i-1,j))/(x(i) - x(i-1)) + U;
          v(z) = (phi(i,j) - phi(i-1,j))/(y(i) - y(i-1));
          UE(z) = sqrt(u(z)^2 + v(z)^2)/U;
          
      
      z = z+1;
end

z=1;
for i = 21:70
    xnew(z) = x(i) - 20;
    z = z+ 1;
end

export = [Nx_2,Ny,xnew,y1(21:70),UE,3000000.00000,0.09920,8.00000,1.14000,0.01000,1.00000];
exp1 = transpose(export);


% figure;
% plot(linspace(20,21,50),cp);
% set(gca,'ydir','reverse'); 
% xlabel('X - axis \rightarrow');
% ylabel('Coefficient of Pressure, C_p \rightarrow');
% title('Coefficient of Pressure along Airfoil Surface, M_\infty = 0.88');
