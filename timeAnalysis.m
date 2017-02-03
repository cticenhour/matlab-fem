% Time analysis for electromagnetic waves solved in frequency domain in FEM MATLAB waveguide code

function [E_time,E_analytic_time] = timeAnalysis(U,node_list,omega,beta,width,init_E,rF_cycles)

t_max = rF_cycles/omega;

t = 0:t_max/100:t_max;

E_time = U*exp(-1i*omega*t);

E_analytic_time = zeros(length(node_list(:,1)),length(t));

for i=1:length(t)    
    E_analytic_time(:,i) = init_E*sin((pi/width)*node_list(:,2)).*cos(omega*t(i)-beta*node_list(:,1));
end

end