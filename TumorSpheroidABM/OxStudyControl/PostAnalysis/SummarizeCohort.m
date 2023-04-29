% a script to summarize cohort data. not sure why I originally created this
% as a function...but it seems more convenient to have this as a script.
% saves mean cell counts, phase counts, ode state variable counts, and all their
% SDs.

clearvars;
cohort_name = "cohort_230124175743017";
addpath("~/Documents/MATLAB/myfunctions/")
load(sprintf("../../data/%s/output.mat",cohort_name),"ids","nsamps_per_condition","cohort_size","lattice_parameters");
nsamps_per_parameter_vector = nsamps_per_condition;
n_conditions = 1;
C = {[]};
vals = {lattice_parameters.values};
%%
for i = numel(ids):-1:1
    S = load(sprintf("../../data/sims/%s/output_final.mat",ids(i)));
    phase_count(:,:,:,i) = reshape(S.tracked.phases,[],2,2);
end

t_abm = S.tracked.t;
t_abm = round(1440*t_abm)/1440; % to make sure that the last time point is actually 3 days (not 3-eps() days)
nt_abm = length(t_abm);

t = [0;10;24;48;72]/24;
nt = length(t);


%% reshape and summarize
phase_count = reshape(phase_count,[nt_abm,2,2,size(ids)]);
phase_count_sampled = interp1(t_abm,phase_count,t);
state_vars = sum(phase_count_sampled,2);
state_vars = reshape(state_vars,nt,2,[],nsamps_per_parameter_vector);
temp_avg = mean(state_vars,ndims(state_vars));
temp_std = std(state_vars,[],ndims(state_vars));

temp_avg = reshape(temp_avg,nt,2,[]);
temp_std = reshape(temp_std,nt,2,[]);

%% create data struct and replace 0 SD with nonzero SD
if isempty(gcp('nocreate'))
    parpool("Threads")
end
r_threshold = 1e-14;
for i = size(temp_avg,3):-1:1
    D(i).A = temp_avg(:,:,i);
    covs = zeros(2,2,nt);
    r = zeros(nt,1);
    parfor ti = 1:nt
        covs(:,:,ti) = cov(squeeze(state_vars(ti,:,i,:))');
        r(ti) = rcond(squeeze(covs(:,:,ti)));
    end

    bad_r = r < r_threshold;
    I = find(bad_r);
    for ii = length(I):-1:1
        if I(ii)==1 && ~bad_r(2) % if it is the first time point and the second is nonzero, use that
            covs(:,:,1) = covs(:,:,2);
            r(1) = r(2);
            bad_r(1) = bad_r(2);
        elseif I(ii)==nt && ~bad_r(nt-1) % similarly for the last time point
            covs(:,:,nt) = covs(:,:,nt-1);
            r(nt) = r(nt-1);
            bad_r(nt) = bad_r(nt-1);
        elseif all(I(ii)~=[1,nt]) % if it's in the middle...
            if ~any(bad_r(I(ii)+[-1,1])) % take the average of the SDs on either side
                covs(:,:,I(ii)) = 0.5*(covs(:,:,I(ii)-1)+covs(:,:,I(ii)+1));
                r(I(ii)) = rcond(squeeze(covs(:,:,I(ii))));
                bad_r(I(ii)) = r(I(ii)) < r_threshold;
                if bad_r(I(ii))
                    covs(:,:,I(ii)) = eye(2);
                    r(I(ii)) = 1;
                    bad_r(I(ii)) = false;
                end
            elseif ~bad_r(I(ii)-1) % if the one prior is nonzero, use that
                covs(:,:,I(ii)) = covs(:,:,I(ii)-1);
                r(I(ii)) = r(I(ii)-1);
                bad_r(I(ii)) = bad_r(I(ii)-1);
            elseif ~bad_r(I(ii)+1) % if the one after is nonzero, use that
                covs(:,:,I(ii)) = covs(:,:,I(ii)+1);
                r(I(ii)) = r(I(ii)+1);
                bad_r(I(ii)) = bad_r(I(ii)+1);
            else % if all else fails, set it to Identity to use the actual error
                covs(:,:,I(ii)) = eye(2);
                r(I(ii)) = 1;
                bad_r(I(ii)) = false;
            end
        else % then it is at the beginning or end with a zero next to it
            covs(:,:,I(ii)) = eye(2);
            r(I(ii)) = 1;
            bad_r(I(ii)) = false;
        end

        if bad_r(I(ii))
            covs(:,:,I(ii)) = eye(2);
            r(I(ii)) = 1;
            bad_r(I(ii)) = false;
        end

    end

    dets = zeros(nt,1);
    invs = zeros(nt,2,2);
    parfor ti = 1:nt
        dets(ti) = det(covs(:,:,ti));
        invs(ti,:,:) = inv(covs(:,:,ti));
        stds(ti,:) = sqrt(diag(covs(:,:,ti)));
    end
    D(i).C = covs;
    D(i).D = dets;
    D(i).Q = invs;
    D(i).S = stds;
end

% for i = size(temp_avg,3):-1:1
%     D(i).A = temp_avg(:,:,i);
%     [I,J] = find(temp_std(:,:,i)==0);
%     for ii = 1:length(I)
%         if I(ii)==1 && temp_std(2,J(ii),i)~=0 % if it is the first time point and the second is nonzero, use that
%             temp_std(1,J(ii),i) = temp_std(2,J(ii),i);
%         elseif I(ii)==nt && temp_std(nt-1,J(ii),i)~=0 % similarly for the last time point
%             temp_std(nt,J(ii),i) = temp_std(nt-1,J(ii),i);
%         elseif all(I(ii)~=[1,nt]) % if it's in the middle...
%             if all(temp_std(I(ii)+[-1,1],J(ii),i)~=0) % take the average of the SDs on either side
%                 temp_std(I(ii),J(ii),i) = mean(temp_std(I(ii)+[-1,1],J(ii),i));
%             elseif any(temp_std(I(ii)+[-1,1],J(ii),i)~=0) % if one is nonzero, use that
%                 temp_std(I(ii),J(ii),i) = max(temp_std(I(ii)+[-1,1],J(ii),i));
%             else % if all else fails, set it to 1 to use the acutal error
%                 temp_std(I(ii),J(ii),i) = 1; 
%             end
%         else % then it is at the beginning or end with a zero next to it
%             temp_std(I(ii),J(ii),i) = 1;
%         end
%     end
%     D(i).S = temp_std(:,:,i);
% end
D = reshape(D,[1,cohort_size]);

save(sprintf("../../data/%s/summary_short.mat",cohort_name),"D","t","C","cohort_size","nsamps_per_parameter_vector","n_conditions","vals","-v7.3")


%% old version
% count = squeeze(sum(phase_count,2));
% 
% 
% %%
% count = reshape(count,[nt,size(ids)]);
% phase_count = reshape(phase_count,[nt,4,size(ids)]);
% ode_state_count = cat(2,sum(sliceof(phase_count,2,1:2),2),sum(sliceof(phase_count,2,3:4),2));
% 
% %%
% average_count = mean(count,ndims(count));
% phase_average = mean(phase_count,ndims(phase_count));
% ode_state_average = mean(ode_state_count,ndims(ode_state_count));
% 
% %%
% count_std = std(count,[],ndims(count));
% phase_std = std(phase_count,[],ndims(phase_count));
% ode_state_std = std(ode_state_count,[],ndims(ode_state_count));

%%
% save(sprintf("../data/%s/summary.mat",cohort_name),"average_count","phase_average","ode_state_average","count_std","phase_std","ode_state_std","-v7.3")
