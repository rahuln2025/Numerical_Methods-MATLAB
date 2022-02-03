mu_min=2.4; mu_max=4; %range of mu values
n_mu=500; %number of mu pixels
n_x=400; %number of x pixels
mu_edges=linspace(mu_min,mu_max,n_mu+1); %edges of mu pixels
mu=(mu_edges(1:n_mu)+mu_edges(2:n_mu+1))/2; %values of mu on which to perform computation
x_edges=linspace(0,1,n_x+1); %edges of x pixels

n_trans=20000; %transient iterations
n_data=10000;  %number of x values per mu value

x_data=zeros(n_data,n_mu); %x-data used to construct figure

x_0=0.5; %initial condition

% WRITE THE COMPUTATIONAL ENGINE OF THE CODE.
% USE THE ALREADY DEFINED PARAMETERS AND VARIABLES: n_mu, mu, x_0, n_trans, n_data.
% YOUR FINAL RESULT WILL BE THE VARIABLE x_data, and this variable will be assessed.
for ite = 1:n_mu
    xold = x_0;
    mu_calc = mu(ite);
    %ite_trans = 1;
    for ite_trans = 1: n_trans
    % ite_trans < n_trans
        xnew = mu_calc*xold*(1 - xold);
        xold = xnew;
        %ite_trans = ite_trans + 1;
    end
    
    %ite_data = 1;
    for ite_data = 1 : n_data
    %while ite_data < n_data
        xnew = mu_calc*xold*(1 - xold);
        x_data(ite_data, ite) = xnew;
        xold = xnew;
        %ite_data = ite_data + 1;
    end
end

%%%%% bin data and plot image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_histogram=zeros(n_x,n_mu); %binned values of x
for i=1:n_mu
x_histogram(:,i)=histcounts(x_data(:,i),x_edges);
x_histogram(:,i)=255*x_histogram(:,i)/max(x_histogram(:,i));
end
colormap(flipud(gray(256))); brighten(-0.8); cmap=colormap;
im=image([mu_edges(1) mu_edges(end)], [x_edges(1) x_edges(end)], x_histogram);
set(gca,'YDir','normal');
xlabel('$\mu$','Interpreter','latex','FontSize',14);
ylabel('$x\;\;$','Interpreter','latex','FontSize',14);
xlim([3 4]);
title('Logistic Map Bifurcation Diagram','Interpreter','latex','FontSize',16)
