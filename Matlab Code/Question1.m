Cn = 0;
Ca = 0;
j=1;
for i=21:70
    dx = x(i+1) - x(i);
    dy = y1(i+1) - y1(i);
    denom = sqrt(dx^2 + dy^2);
    nx = -dy/denom;
    ny = dx/denom;
    
    Cn = Cn + (-cp(i-20)*ny*denom);
    Ca = Ca + (-cp(i-20)*nx*denom);
end

Cdp90 = 2*(Cn*sind(0) + Ca*cosd(0))
Clp = Cn*cosd(0) - Ca*sind(0);