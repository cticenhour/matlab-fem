% Function to perform basic electromagnetic wave analysis for FEM MATLAB code. OUTPUTS: [k_axis,spectrum,k_parallel,R,k_calc,power]

function [k_axis,spectrum,R,R_theory,k_calc,power] = FFTanalysis(node_list,U,len,width,num_pts,method,windowed)

location = width/2;

xy = [linspace(0,len,num_pts)', ones(num_pts,1).*(location)];
z = griddata(node_list(:,1),node_list(:,2),U,xy(:,1),xy(:,2),method);

if windowed == 1
    window_width = num_pts*0.9;
    w = hann(window_width);
    w = [zeros((num_pts - window_width)/2,1); w; zeros((num_pts - window_width)/2,1)];
    z = w.*z;
end

k_axis = 2.*pi./mean(diff(linspace(0,len,num_pts))).*[-num_pts/2:(num_pts/2-1)]./num_pts;

spectrum = fftshift(abs(fft(z')/len));


halfway = round(length(spectrum)/2);
[fwd_mag,fwd_index] = max(spectrum);
[~,rev_index] = min(sqrt(k_axis(fwd_index)^2 - k_axis(1:halfway).^2));

rev_mag = spectrum(rev_index);      % Potential issue here. Neglects when reflection magnitude is drowned out by broad forward wave peak.

% Calculate reflection coefficient
R = rev_mag/fwd_mag;

% Calculate k-components of forward traveling wave 
k_parallel = k_axis(fwd_index);
k_perp = pi/width;
k_calc = sqrt(k_perp^2 + k_parallel^2);

% Calculate wave incident angle on surface
angle = atan(k_perp/k_parallel);

% Calculate theoretical reflection coefficient
R_theory = abs((cos(angle) - 1)/(cos(angle) + 1));
%R_theory = abs((cos(angle) + 0.5*sin(angle)^2 - 1)/(cos(angle) - 0.5*sin(angle)^2 + 1));

power = abs(spectrum).^2;

end