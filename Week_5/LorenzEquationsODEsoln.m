sigma=10; beta=8/3; r=28;
x0=1; y0=1; z0=1; tspan=[0 100];
ntrans=20;
options = odeset('RelTol',1.e-6);
[t,xyz]=ode45(@(t, xyz) lorenz_eqs(xyz,sigma,beta,r), tspan, [x0, y0, z0], options);
x=xyz(ntrans:end,1); y=xyz(ntrans:end,2); z=xyz(ntrans:end,3);
plot3(x,y,z); 
xlabel('$x$','Interpreter','latex','FontSize',14 );
ylabel('$y$','Interpreter','latex','FontSize',14 );
zlabel('$z$','Interpreter','latex','FontSize',14 );
title('Lorenz Equations','Interpreter','latex','FontSize',16);

function dxyzdt = lorenz_eqs(xyz,sigma,beta,r)
x=xyz(1); y=xyz(2); z=xyz(3);
dxyzdt=[sigma*(y - x); x*(r - z) - y; x*y - beta*z];
end
