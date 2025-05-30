%%%%%%%%%%%%%%%%%%%% DEPRECATED %%%%%%%%%%%%%%%%%%%%%%%%
% This script sets up and calls the profile likelihood method for
% identifiability for the ODE pars from the data. it then plots the
% combinations

clearvars;

addpath("~/Documents/MATLAB/myfunctions/") % replace with path (rel or abs) to myfunctions
addpath("../ODEFitting/")

% folder to store plots and text files

% load data and compute the standard error
ODE_fit = load("../ODEFitting/ODEFittoData.mat");

%%
threshold = chi2inv(0.95,3); % compute threshold value for the parameter confidence intervals

% specify parameter ranges
para_ranges = [0,1e1;     % lambda
               0,1e1;  % alpha
               0,2e3];      % K

lb = [0;0;0];
ub = [1e1;1e1;2e3];
opts = optimset('Display','off','TolFun',1e-12,'TolX',1e-12);

% specify parameter names
para_names = {'\lambda', '\alpha', 'K'};
para_names_file_save = {'lambda', 'alpha', 'K'};

N = size(para_ranges,1); % number of parameters

%% compute profile likelihoods
out = profileLikelihood(ODE_fit.pstar,ODE_fit.tt,ODE_fit.data,ODE_fit.data_std,para_ranges,lb,ub,opts,threshold,true);

%% save the output
% save("ProfileLikelihoods_Data.mat");
% ProfileLikelihoods_DataRestricted.mat used better bounds for the
% parameters in fitting. Use that.

%%  plot parameter combinations
figure;
ci=0;
for i = 1:N-1
    for j = i+1:N
        ci = ci+1;
        subplot(2,N,r2c(2,N,[1,ci]))
        plot(out{i}(i,:),out{i}(j,:))
        xlabel(para_names{i});ylabel(para_names{j})
        subplot(2,3,r2c(2,3,[2,ci]))
        plot(out{j}(j,:),out{j}(i,:))
        xlabel(para_names{j});ylabel(para_names{i})
    end
end

%%  plot parameter combinations
figure;
ci=0;
for i = 1:N-1
    for j = i+1:N
        ci = ci+1;
        subplot(1,N,r2c(1,N,[1,ci])); hold on
        plot(out{i}(i,:),out{i}(j,:)) % parameter j values as i was profiled
        xlabel(para_names{i});ylabel(para_names{j})
        plot(out{j}(i,:),out{j}(j,:)) % parameter i values as j was profiled (plotted as para i vs para j so it shares the same axes)
        xlabel(para_names{j});ylabel(para_names{i})
    end
end

%%

rmpath("../ODEFitting/")
