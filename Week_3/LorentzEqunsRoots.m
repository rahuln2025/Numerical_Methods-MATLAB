r=28; sigma=10; b=8/3; 
RelTol=1.e-06; AbsTol=1.e-09;
for nroot=1:3
    if nroot==1, x=1; y=1; z=1; end
    if nroot==2, x=10; y=10; z=10; end
    if nroot==3, x=-10; y=-10; z=-10; end
    error=Inf;
    while error > max(RelTol*max(abs([x,y,z])),AbsTol)
        J= [-sigma, sigma, 0; r - z, -1, -x; y, x, -b];
        rhs = -[sigma*(y - x); x*(r - z) - y; x*y - b*z]; 
        delta_xyz=J\rhs;
        x = x + delta_xyz(1);
        y = y + delta_xyz(2);
        z = z + delta_xyz(3);
        error=max(abs(delta_xyz));
    end
    xroot(nroot)=x; yroot(nroot)=y; zroot(nroot)=z;
end
roots=[xroot;yroot;zroot];
fprintf('steady-state solution:\n')
fprintf('(x, y, z) = (%2.0f,%2.0f,%2.0f)  \n', roots(:,1));
fprintf('(x, y, z) = (%7.5f,%7.5f,%3.0f)  \n', roots(:,2));
fprintf('(x, y, z) = (%7.5f,%7.5f,%3.0f)  \n', roots(:,3));
