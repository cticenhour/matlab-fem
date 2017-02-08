% MATLAB function of "Analytic" series solution of U_xx + U_yy = 1 on a
% unit square with U(x,y) = 0 at the boundaries
% Author: Casey Tyler Icenhour
% Organization: NCSU/ORNL
% September 26, 2016

function U_analytic = analyticFXN_poisson2D(p,series_index_max)

% generate ID list of nodes from length of p X-Y array
nodes = 1:1:length(p);

% Stopping the infinite series at a given index in x and y components to 
% approximate the solution - hence "analytic"
%series_index_max = 30;

U_analytic = zeros(length(nodes),1);

k = 1;
l = 1;
m = 1;

while m < length(nodes)
    x = p(nodes(m),1);
    y = p(nodes(m),2);
    while k < series_index_max
        while l < series_index_max
            U_analytic(m) = U_analytic(m) + 4*(1-cos(pi*k))*(1-cos(pi*l))*sin(k*pi*x)*sin(l*pi*y)/(pi^4*k*l*(k^2+l^2));
            l=l+1;
        end
        l = 1;
        k=k+1;
    end
    k = 1;
    m=m+1;
end