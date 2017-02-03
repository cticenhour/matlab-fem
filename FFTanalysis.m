% Function to perform basic electromagnetic wave analysis for FEM MATLAB code. OUTPUTS: [k_axis,spectrum,k_parallel,R,k_calc,power]

function [k_axis,spectrum,k_parallel,R,k_calc,power] = FFTanalysis(k,node_list,U,len,width,num_pts,method,windowed)

location = width/2;

xy = [linspace(0,len,num_pts)', ones(num_pts,1).*(location)];
z = griddata(node_list(:,1),node_list(:,2),U,xy(:,1),xy(:,2),method);

if windowed == 1
    window_width = num_pts*0.9;
    w = hann(window_width);
    w = [zeros((num_pts - window_width)/2,1); w; zeros((num_pts - window_width)/2,1)];
    z = w.*z;
end

k_axis = (-num_pts/2 + 1:num_pts/2)*k;

spectrum = fftshift(abs(fft(z)));


halfway = round(length(spectrum)/2);
[fwd_mag,fwd_index] = max(spectrum(halfway:length(spectrum)));
rev_mag = max(spectrum(1:halfway));
% Calculate reflection coefficient
R = rev_mag/fwd_mag;
%
k_parallel = k_axis(fwd_index);
k_perp = pi/width;
k_calc = sqrt(k_perp^2 + k_parallel^2);

power = abs(spectrum).^2;

end