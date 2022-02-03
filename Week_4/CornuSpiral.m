c = @(x) cos(0.5*pi*(x.^2)); % assign the integrand for C(t)
s = @(x) sin(0.5*pi*(x.^2)) % assign the integrand for S(t)

tmin=-8; tmax=8; nt=2000;
t=linspace(tmin,tmax,nt);
C=zeros(nt,1); S=zeros(nt,1);
for i=1:nt
    C(i)=integral(@(x) c(x), 0, t(i)); % compute C(i) using integral.m and the integrand c(x) defined on top
    S(i)=integral(@(x) s(x), 0, t(i)); % compute S(i) using integral.m and the integrand s(x) defined on top
end
plot(S,C)
xlabel('$S(t)$', 'Interpreter', 'latex', 'FontSize',14);
ylabel('$C(t)$', 'Interpreter', 'latex', 'FontSize',14);
title('Cornu spiral', 'Interpreter', 'latex', 'FontSize',16);