% Compute the Feigenbaum delta
% Store approximate values in the row vector delta for assessment, where length(delta)= num_doublings and 
% delta(2:num_doublings) are computed from the algorithm described in Lectures 21-23.
num_doublings=11; delta=zeros(1,num_doublings+1); delta(1)=5;
% Write your code here

%compute initial guess for m value to start with Newton's method
% m  stores the guesses of m's
m = zeros(1, num_doublings);
m(1) = 2; %m0
m(2) = 1 + sqrt(5); %m1
% m_values stores the correct m values
m_values = zeros(1, num_doublings);
m_values(1) = 2;%m0
m_values(2) = 1 + sqrt(5); %m1

for n = 2:num_doublings    
    %initial guess of m(n) for Newton's method
    m(n+1) = m(n) + (m(n) - m(n-1))/delta(n-1) ;
    
    for nit = 1:40
        %iterate Newton's method 40 times (should be enough)
        X = 0.5;
        X_prime = 0;
        N = 2^(n);
        for i = 1:N
            %iterate till the period is complete i.e. x0 to xn
            %find Xn and X'n
            X_temp = m(n+1) * X * (1 - X);
            X_prime_temp = X * (1 - X) + (m(n+1) * X_prime * (1 - 2*X));
            
            X = X_temp;
            X_prime = X_prime_temp;
            
        end
        %iterate Newton's method to find m (growth rate)
        m(n+1) = m(n+1) - (X - 0.5) / X_prime;    
    end
    %save values of m and find delta from previous values
    m_values(n+1) = m(n+1);
    delta(n) = (m_values(n) - m_values(n-1)) / (m_values(n+1) - m_values(n));
end


% Output your results
fprintf('n        delta(n)\n');
for n=1:num_doublings
    fprintf('%2g %18.15f\n',n,delta(n));
end