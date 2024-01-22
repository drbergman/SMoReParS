clearvars;

% this script will test out how to do sensitivity on the ABM using the SM

%% Program to run

addpath("~/Documents/MATLAB/myfunctions/")
addpath("../ODEFitting/")
addpath("../ProfileLikelihood/")
addpath("../../../SensitivityFns/")
addpath("../../../ProfileLikelihoodFns/")

use_profiles = true;
%
% This algorithm is an adaptation of the method of Sensitivity Analysis
% called the Morris method.
%
% Sensitivity analysis is used to estimate the influence of uncertainty
% factors on the output of a function.
% The Morris method is sometimes referenced to as a qualitative method : it
% gives rough estimations with a limited number of calculations.
% The Morris method can be used to simplify a function, as a first step. It
% can identify the factors with a low influence which can be fixed.
% For further information :
% Saltelli, A., Tarantola, S., Campolongo, F., and Ratto, M. (2004). Sensitivity Analysis in Practice - A Guide to Assessing Scientific Models. Wiley.
% 
% This algorithm reduces the risk to underestimate and fix non-negligible
% factors. It is presented in:
% Henri Sohier, Helene Piet-Lahanier, Jean-Loup Farges, Analysis and optimization of an air-launch-to-orbit separation, Acta Astronautica, Volume 108, March泡pril 2015, Pages 18-29, ISSN 0094-5765, http://dx.doi.org/10.1016/j.actaastro.2014.11.043.
% (http://www.sciencedirect.com/science/article/pii/S0094576514004974)
%
% This program is divided in 6 parts:
% 1) Clearing the memory
% 2) Parameters : Please fill in
% 3) Initialization of the variables
% 4) Loop
% 5) Output text
% 6) Figure
%
% Please fill in the second part to apply the algorithm to your function.
% Do not change the parameters to see the results with the "modified Sobol
% test function".
%
% This program outputs a figure as well as a short summary in the console.
% Consider fixing the factors which appear as negligible on the left of the
% figure (not necessary all the factors under the limit).
%% 2) Parameters : Please fill in
% Maximum number of simulation runs :
% Large number = better estimation of the influence of the factors
% Recommended value : (number of factors + 1) * 10
% The algorithm will maybe exceed this value if it is considered necessary
options.limit_factor = 0.5; % how to set the limit for separating low and high impact factors
options.initialization_factor = 0.8; % how to determine when it has been sufficiently initialized
nsim_max = 1e4;
nsamps = 100; % number of points to sample in LHS for ODE pars
% Function studied :
% Replace test_function by the name of your function. It must be a 
% function with one multidimensional input x. x must represent the values 
% of the uncertainty factors in the quantiles hyperspace (the i-th 
% coordinate of x is not the actual value of the i-th factor, but the 
% corresponding value of the cumulative distribution function of the i-th 
% factor). To adapt your function, first calculate the actual values of 
% the factors by applying the inverse of their cumulative distribution 
% function to each coordinate of x; Matlab includes such inverses: 
% mathworks.com/help/stats/icdf.html ).

% par_names = ["carrying_capacity";"occmax_2d";"move_rate_microns";"g1_to_s";"s_to_g2";"g2_to_m";"m_to_g1"];
files.profiles = "../ProfileLikelihood/data/Profiles_SMFromABM_New_clean.mat";
D = makeABMParameterDistributionsDictionary();

%% create bounding hypersurfaces


% Number of factors of uncertainty of the function studied :
% num_factors=numel(par_names); 

% assert(num_factors==numel(par_names)) % make sure that there is a value for each of the parameters to be varied
% assert(D.numEntries==numel(par_names)) % make sure each parameter has an associated distribution
T = dictionary("occmax_2d", @(x) min(7,floor(x)));
sm_functional = @(p) sum(computeTimeSeries(p,[],[],false,3));
% studied_function = setupSampleFromSMFunction(files,sm_functional,par_names=par_names, T = T,D=D,nsamps=nsamps,use_profiles=use_profiles);
[studied_function,num_factors] = setupSampleFromSMFunction(files,sm_functional,T=T,D=D,nsamps=nsamps,use_profiles=use_profiles);
[mu_star,sigma,~] = morris_simple(studied_function,num_factors,15,sort_output=false);

%% 3) Initialization of the variables
% table_outputs = []; % All the outputs of the simulations runs. One line = results around one point of the factors hyperspace. First column = output at a sampled point, second column = output after varying the first factor, etc...
% table_ee = []; % All the elementary effects. One line = elementary effects at one point of the factors hyperspace. First column = elementary effect of the first factor, etc... See the relation between the outputs of the simulation run and the elementary effects.
% factors_over = []; % Indexes of the factors over the limit (important factors). The elementary effects of these factors will not be calculated at the next step.
% table_factors_over = []; % Factors over the limit at the different steps. n-th line = factors with no elementary effect at the n-th step = factors over the limit at the (n-1)-th step.
% points=[]; % Sampled points of the factors hyperspace where the elementary effects are calculated.
% n=1; % Current step.
% nsim = num_factors+1; % Number of simulation runs after the next step.
% initialization = 0; % Boolean, the calculations will foccus on the factors under the limit when it will equal 1.
% convergence = 0; % Boolean, equals 1 when the factors over the limit have not changed over the last steps.
% %% 4) Loop
% while (nsim <= nsim_max) || (n<=size(table_ee,1)) || ~convergence 
% % Continues if at least one of the conditions is true
% % Condition 1 : The calculations can continue if the number of simulation runs after the next step (nsim) is not larger than the maximum number of simulation runs
% % Condition 2 : The calculations can continue if the algorithm returned to a previous step and did not finish to complete the table of elementary effects
% % Condition 3 : The calculations can continue if the set of factors over the limit has changed over the last steps
%     % Application of the algorithm at the current step:
%     [ table_outputs, table_ee, factors_over, n, points, initialization ] = moat_loop_function( num_factors, studied_function, table_outputs, table_ee, factors_over, n, points, initialization );
%     % Updating values
%     nsim = sum(table_ee(:)~=-1)+size(table_ee,1)+length(factors_over)+1;
%     table_factors_over(n,1:length(factors_over)) = factors_over;
% 
%     % Comparison of the factors over the limit over the three last steps
%     if n>=4 % After one application of moat_loop_function, n=2. Thus, three steps are completed when n>=4.
% 
%         % Last step
%         factors_n = table_factors_over(n,:); % n-th line of table_factors_over
%         to_delete = factors_n==0; % Depending on the number of factors over the limit at the other steps, there may be zeros in factors_n
%         factors_n(to_delete) = []; % The zeros are deleted (they are not factors indexes).
%         factors_n = sort(factors_n); % The factors are sorted in ascending order (to compare the sets of factors regardless of their values/orders)
%         % Last step - 1
%         factors_n_1 = table_factors_over(n-1,:);
%         to_delete = factors_n_1==0;
%         factors_n_1(to_delete) = [];
%         factors_n_1 = sort(factors_n_1);
%         % Last step - 2
%         factors_n_2 = table_factors_over(n-2,:);
%         to_delete = factors_n_2==0;
%         factors_n_2(to_delete) = [];
%         factors_n_2 = sort(factors_n_2);
%         convergence = isequal(factors_n,factors_n_1) & isequal(factors_n_1,factors_n_2); % Equals 1 when the factors over the limit have not changed over the three last steps
%     else
%         convergence = 0;
%     end
% end
% %% 5) Output text
% disp('*******************************');
% disp('* SUMMARY OF THE CALCULATIONS *');
% disp('*******************************');
% disp(['Number of factors : ' num2str(num_factors)]);
% disp(['Chosen number of simulation runs : ' num2str(nsim_max)]);
% disp(['Actual number of simulation runs : ' num2str(sum(table_ee(:)~=-1)+size(table_ee,1))]);
% disp(['Number of points tested in the hyperspace : ' num2str(size(table_ee,1))]);
% disp(['Number of points normally tested with the same number of simulation runs  : ' num2str(floor((sum(table_ee(:)~=-1)+size(table_ee,1))/num_factors))]);
% disp('*******************************');
% %% 6) Figure
% max_ee = max(table_ee,[],1); % Maxima of the elementary effects of the factors.
% [sort_A, sort_B] = sort(max_ee,'ascend'); % Ordering the maxima.
% sort_Ab = [sort_A(1) sort_A(1:num_factors-1)]; % New table where the first elementary effect is repeted twice.
% difference = sort_A-sort_Ab; % Variation between the successive maxima.
% dmax = max(difference); % Largest variation.
% sep = dmax*options.limit_factor; % Amplitude of the variation which can be considered as the limit.
% sep_indic = find(difference>=sep, 1 ); % Index of the first factor after the limit
% % Opening the figure
% hfig = figure;
% hold on;
% for j=1:num_factors % For each factor
%    index_fac = sort_B(j);
%    plot(j,table_ee(table_ee(:,index_fac)~=-1,index_fac),'*','LineWidth',2);
% %    for k=1:size(table_ee,1) % For each elementary effect
% %        y = table_ee(k,index_fac); % Value of the elementary effect.
% %        if y~=-1 % -1 represents a value which has not been calculated
% %            plot(j,y,'*','linewidth',2); % Elementary effects.
% %        end
% %    end
% end
% plot([sep_indic-0.5 sep_indic-0.5], [0 1.1*max(table_ee(:))], 'k', 'linewidth', 2); % Limit.
% hold off;
% % Factors indexes on the x-axis
% labels = cell(num_factors,1);
% for l=1:num_factors
%    labels{l} = num2str(sort_B(l));
% end
% set(gca, 'XTick',1:num_factors, 'XTickLabel', labels, 'FontSize',12)
% axis([0 num_factors+1 0 1.1*max(table_ee(:))]) % Limits of the axes.
% % Labels
% xlabel('Factors ordered by ascending maximum','FontSize',12)
% ylabel('Elementary effects','FontSize',12)


rmpath("../ODEFitting/")
rmpath("../ProfileLikelihood/")
rmpath("../../../SensitivityFns/")
rmpath("../../../ProfileLikelihoodFns/")