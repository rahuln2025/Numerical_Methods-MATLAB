num_roots=5; num_functions=6;
%initial guess for roots (from Wolfram MathWorld)
zeros_guess=[2.4,3.8,5.1,6,7.5,8.7;...
   5.5,7,8.4,9.7,11,12;...
   8.6 10,11.6,13,14,16;...
   11.8,13,15,16,18,19;...
   15,16.4,18,19.4,21,22]; 
%Compute first five roots of first six Bessel functions
%Put in variable bzeros with size(bzeros) = [5, 6]
%bzeros = zeros(5,6);
%x = linspace(0, 25, 100);
%theta = linspace(0, pi, 100);
for n = 1:num_functions %functions
    fun = @(theta,x,n) cos(x.*sin(theta) - (n-1)*(theta));
    %J = (1/pi)*(integral(@(theta) fun(theta, x, n), 0, pi, 'ArrayValued', 1));
    %plot(x, J)
    fxn = @(x) (1/pi)*(integral(@(theta) fun(theta, x, n), 0, pi, 'ArrayValued', 1)) ;
    
    for k = 1:5 %roots
        x0 = zeros_guess(k,n);
        bzeros(k,n) = fzero(@(x) fxn(x) ,x0);
    end
end

%print table
fprintf('k     J0(x)     J1(x)     J2(x)     J3(x)     J4(x)     J5(x)\n')
for k=1:num_roots
    fprintf('%i',k)
    for n=0:num_functions-1
        fprintf('%10.4f',bzeros(k,n+1));
    end
    fprintf('\n');
end
