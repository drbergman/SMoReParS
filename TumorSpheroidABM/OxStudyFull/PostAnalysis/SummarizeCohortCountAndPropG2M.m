% a script to summarize cohort data just with count and proportion in g2/m.

clearvars;
cohort_name = "cohort_2303301105";
addpath("~/Documents/MATLAB/myfunctions/")

%% load cohort data
cohort = load(sprintf("../../data/%s/output.mat",cohort_name));
nsamps_per_parameter_vector = cohort.nsamps_per_condition;
vals = {cohort.lattice_parameters.values};

%% get conditions
n_conditions = 3; % [control;0.75;7.55] uM of Oxaliplatin
LP_paths = arrayify(cohort.lattice_parameters,"path");
condition_dim = find(all(LP_paths == ["chemo_pars","concentration"],2));
C = mat2cell(cohort.lattice_parameters(condition_dim).values,[1 1 1]);

%% load tracked data and concatenate
for i = numel(cohort.ids):-1:1
    S = load(sprintf("../../data/sims/%s/output_final.mat",cohort.ids(i)));
    if size(S.tracked.phases,2)>4
        error("Not sure how we will count the arrested compartment in this.")
    end
    phase_count(:,:,:,i) = reshape(S.tracked.phases,[],2,2);
end

t_abm = S.tracked.t;
t_abm = round(1440*t_abm)/1440; % to make sure that the last time point is actually 3 days (not 3-eps() days)
nt_abm = length(t_abm);

t = [10;24;36;48;72]/24; % do not bother comparing the first time point since none of the varied ABM parameters affect the initial conditions
nt = length(t);

%% reshape and summarize
phase_count = reshape(phase_count,[nt_abm,2,2,size(cohort.ids)]);
phase_count_sampled = interp1(t_abm,phase_count,t); % All this does is remove the t=0 point as the ABM data was saved just at the sample data points
count = sum(phase_count_sampled,2:3);
count = reshape(count,nt,[],nsamps_per_parameter_vector);
count_avg = mean(count,ndims(count));
count_std = std(count,[],ndims(count));

state_vars = sum(phase_count_sampled,2);
state_vars = reshape(state_vars,nt,2,[],nsamps_per_parameter_vector);
state2_prop = squeeze(state_vars(:,2,:,:))./count;
state2_prop_avg = mean(state2_prop,ndims(state2_prop));
state2_prop_std = std(state2_prop,[],ndims(state2_prop));
% 
% all_avg = mean(state_vars,ndims(state_vars));
% all_std = std(state_vars,[],ndims(state_vars));
% 
% all_avg = cat(3,count_avg,state2_prop_avg);
% all_std = cat(3,count_std,state2_prop_std);
% all_std = reshape(all_std,nt,2,[]);

%% create data struct and replace 0 SD with nonzero SD
% if isempty(gcp('nocreate'))
%     parpool("Threads")
% end
for i = size(count_avg,2):-1:1
    D(i).A = [count_avg(:,i),state2_prop_avg(:,i)];
    temp_std = [count_std(:,i),state2_prop_std(:,i)];
    [I,J] = find(temp_std==0);
    for ii = 1:length(I)
        if I(ii)==1 && temp_std(2,J(ii))~=0 % if it is the first time point and the second is nonzero, use that
            temp_std(1,J(ii)) = temp_std(2,J(ii));
        elseif I(ii)==nt && temp_std(nt-1,J(ii))~=0 % similarly for the last time point
            temp_std(nt,J(ii)) = temp_std(nt-1,J(ii));
        elseif all(I(ii)~=[1,nt]) % if it's in the middle...
            if all(temp_std(I(ii)+[-1,1],J(ii))~=0) % take the average of the SDs on either side
                temp_std(I(ii),J(ii)) = mean(temp_std(I(ii)+[-1,1],J(ii)));
            elseif any(temp_std(I(ii)+[-1,1],J(ii))~=0) % if one is nonzero, use that
                temp_std(I(ii),J(ii)) = max(temp_std(I(ii)+[-1,1],J(ii)));
            else % if all else fails, set it to 1 to use the acutal error
                temp_std(I(ii),J(ii)) = 1; 
            end
        else % then it is at the beginning or end with a zero next to it
            temp_std(I(ii),J(ii)) = 1;
        end
    end
    D(i).S = temp_std;
end

%% reshape to match [conditions,cohort_size]
D = reshape(D,cohort.cohort_size);

cohort_size = cohort.cohort_size(setdiff(1:ndims(D),condition_dim));

D = permute(D,[condition_dim,cohort_size]);

n_time_series = size(D(1).A,2);
save(sprintf("../../data/%s/summary_new.mat",cohort_name),"D","t","C","cohort_size","nsamps_per_parameter_vector","n_conditions","vals","n_time_series","condition_dim","-v7.3")


% %% old code
% count = reshape(count,[nt,size(cohort.ids)]);
% count_average = mean(count,ndims(count));
% count_std = std(count,[],ndims(count));
% 
% phase_count = reshape(phase_count,[nt,2,2,size(cohort.ids)]); % [time, [phase], sample] where [phase] has the 4 phases in this format [G1,G2;S,M] so summing along columns gives the two halves G1/S and G2/M
% ode_state_count = squeeze(sum(phase_count,2));
% state2_prop = squeeze(sliceof(ode_state_count,2,2))./count;
% state2_prop_mean = mean(state2_prop,ndims(state2_prop));
% state2_prop_std = std(state2_prop,[],ndims(state2_prop));
% 
% %%
% filename = sprintf("../../data/%s/summary.mat",cohort_name);
% if exist(filename,"file")
%     save(filename,"count_*","state2_prop_*","-append")
% else
%     save(filename,"count_*","state2_prop_*","-v7.3")
% end
