clear all

c = 3e8;                % m/s

f = 20e6;               % Hz (1/s)
omega = 2*pi*f;         % rad/s

k = omega/c;            % 1/m

a = 8;              % m - Waveguide width

k_c = pi/a;         % 1/m - mode number = 1

beta = sqrt(k^2 - k_c^2);

wavelength = 2*pi/beta; % m

four_wavelengths = 4*wavelength
half_wavelength = wavelength/2