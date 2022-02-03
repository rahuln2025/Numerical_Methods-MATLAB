close all, clear all   clc
%define functions and derivative
f = @(z) z.^4 - 1; fp = @(z) 4*z.^3; %note the elementwise operators used becuase z is an array

%define roots
r1 = 1;
r2 = 1i*1;
r3 = -1;
r4 = 1i*(-1);

%set domain and grid
nx = 10000; ny = 10000;
xmin = -2; xmax = 2; ymin = -2; ymax = 2;
x = linspace(xmin, xmax, nx); y = linspace(ymin, ymax, ny);
[X, Y] = meshgrid(x, y);
Z = X + 1i*Y; %define each point in the domain in complex form

%root finding by Newton's method
nit = 40;
for n = 1:nit
    Z = Z - f(Z) ./ fp(Z);
end

%identify to which root each point converges
eps = 0.001;
Z1 = abs(Z - r1) < eps;
Z2 = abs(Z - r2) < eps;
Z3 = abs(Z - r3) < eps;
Z4 = abs(Z - r4) < eps;
Z5 = ~ (Z1 + Z2 + Z3 + Z4);

%plot figure
map = [ 1 0 0; 0 1 0; 0 0 1; 0 1 1; 0 0 0]; colormap(map); %[red, green, blue, ran color black]
Z_plot = (Z1 + 2*Z2 + 3*Z3 + 4*Z4 + 5*Z5);
image([xmin xmax], [ymin ymax], Z_plot);
set(gca, 'YDir', 'normal');

%beautify plot
axis equal; axis tight;

