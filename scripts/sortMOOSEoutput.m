% sortMOOSEoutput.m 
% Script to sort MOOSE output to fit MATLAB test code node organization
%
% Casey Icenhour
% November 29, 2016
%
% NOTE: Assumes that MATLAB FEM code run data is already present in
% workspace!
% NOTE2: Also assumes that MOOSE code and MATLAB code are using same mesh!

name_query = 'What is the name of the MOOSE output CSV file? ';
filename = input(name_query,'s');
TF = length(strfind(filename, '.csv'));
TF2 = length(strfind(filename, '.CSV'));

if TF==0 && TF2==0
    filename = strcat(filename,'.csv');
end

% MOOSE data import
data = csvread(filename,1,0);

% required arrays initialization
U_MOOSE = zeros(length(data(:,1)),1);

node_list_MOOSE_temp = data(:,4:5);

dist = zeros(length(U_MOOSE),1);

% Sorting loops based on distance functions
for i=1:length(node_list_MOOSE_temp)
    dist = sqrt((node_list(:,1) - node_list_MOOSE_temp(i,1)).^2 + ...
        (node_list(:,2) - node_list_MOOSE_temp(i,2)).^2);
    [min_dist,min_index]=min(dist);
    U_MOOSE(min_index) = data(i,1);
end

