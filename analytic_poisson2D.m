% "Analytic" series solution of U_xx + U_yy = 1 on a
% unit square with U(x,y) = 0 at the boundaries
% Author: Casey Tyler Icenhour
% Organization: NCSU/ORNL
% September 23, 2016

x_min = 0;
y_min = 0;
x_max = 1;
y_max = 1;

interval = 0.05;

x = x_min:interval:x_max;
y = x_min:interval:y_max;

% Stopping the infinite series at a given index in x and y components to 
% approximate the solution - hence "analytic"
series_index_max = 30;

U_analytic = zeros(length(x),length(y));

k = 1;
l = 1;
p = 1;
q = 1;

while k < length(x)
    while l < length(y)
      while p < series_index_max
          while q < series_index_max
              U_analytic(k,l) = U_analytic(k,l) + 4*(1-cos(pi*p))*(1-cos(pi*q))*sin(p*pi*x(k))*sin(q*pi*y(l))/(pi^4*p*q*(p^2+q^2));
              q=q+1;
          end
          q = 1;
          p=p+1;
      end
      p = 1;
      l=l+1;
    end
    l = 1;
    k=k+1;
end