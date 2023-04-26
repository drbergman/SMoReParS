function [time_points,cell_counts,cell_counts_mean,cell_counts_std,cell_counts_se] = data_mean_std()

%% This is the function that provides all the data that the profile 
%% likelihood needs. It is called on the main script and in any others that 
%% require the mean and std data. All data should be provided to the other 
%% files using this file as it makes it easier to run different datasets 
%% this way.

%% Kerri's data
%%{ 
% data time points
time_points = 0:6/24:1794/24;

% ABM parameters
p_div = "0.245";
div_lim = "8";

% Data for Total Cell Number
cell_counts(:,1) = importdata("AllDataCell/AA" + p_div + "AB" + div_lim + "/NumberofCells_CART" + p_div + "_Trial1_t300.txt");    
cell_counts(:,2) = importdata("AllDataCell/AA" + p_div + "AB" + div_lim + "/NumberofCells_CART" + p_div + "_Trial2_t300.txt");
cell_counts(:,3) = importdata("AllDataCell/AA" + p_div + "AB" + div_lim + "/NumberofCells_CART" + p_div + "_Trial3_t300.txt");
cell_counts(:,4) = importdata("AllDataCell/AA" + p_div + "AB" + div_lim + "/NumberofCells_CART" + p_div + "_Trial4_t300.txt");
cell_counts(:,5) = importdata("AllDataCell/AA" + p_div + "AB" + div_lim + "/NumberofCells_CART" + p_div + "_Trial5_t300.txt");
cell_counts(:,6) = importdata("AllDataCell/AA" + p_div + "AB" + div_lim + "/NumberofCells_CART" + p_div + "_Trial6_t300.txt");

cell_counts_mean = mean(cell_counts,2);
cell_counts_std  = std(cell_counts,[],2);
cell_counts_se   = cell_counts_std/sqrt(6);

%}

%% Experimental breast cancer data
%{
cancer_type = "breast";

data = readmatrix("experimental_data/"  + cancer_type + "_cancer_dataset.csv");
std_data = readmatrix("experimental_data/" + cancer_type + "_cancer_std.csv");

time_points = data(:,1);

cell_counts_mean = data(:,2)';
cell_counts_std = std_data(:,2)' - cell_counts_mean;

%}

%{

data = readmatrix("s_shaped_data/data_3/s_shaped_data_mean.csv");
std_data = readmatrix("s_shaped_data/data_3/s_shaped_data_std.csv");

time_points = data(1:end,1);

cell_counts_mean = data(1:end,2)';
cell_counts_std = std_data(1:end,2)' - cell_counts_mean;

%}