nu = U/(3*10^6);
a = [0.25;0.50;0.75];
for j = 1:3
    for i=1:45
        ynew(i,j) = eta(i,j)*sqrt(a(j)/(3*10^6));
    end
end

figure;
plot(U25,ynew(1:39,1))
title('Boundary Layer Profile');
xlabel('u \rightarrow');
ylabel('y-coordinate \rightarrow');

figure
plot(U50,ynew(1:43,2))
title('Boundary Layer Profile');
xlabel('u \rightarrow');
ylabel('y-coordinate \rightarrow');

figure
plot(U75,ynew(1:45,3))
title('Boundary Layer Profile');
xlabel('u \rightarrow');
ylabel('y-coordinate \rightarrow');