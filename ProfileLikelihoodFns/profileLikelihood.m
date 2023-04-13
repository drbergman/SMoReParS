function out = profileLikelihood(pbest,t,D,C,objfn_constants,profile_params,save_all_pars)

% profiles each parameter in pbest. pbest is the best fit at the point. the
% time series to compare against is given in tt, data, and data_std.
% para_ranges sets bounds for profiling. lb and ub set bounds for
% optimization while another parameter is being profiled. opts carries
% fmincon options. threshold sets a cutoff for when profiling can be
% stopped. save_all_pars can be passed in as true to save all parameter
% values, not just the profiled parameter value

% In each direction of pbest(i) (left and right), the ith parameter changes
% at 1% of pbest(i) for 10 steps. If it is still within para_ranges(i,:)
% and the value of the output is within F(pbest)+threshold, then continue
% in that direction until one of those conditions fails.

% will not assume that the pbest came from the same objective function as
% we are using here

m = numel(D); % number of conditions used

% F = @(p) arrayfun(@(j) rawError([p(1:5);objfn_constants.hill_coefficient;p(6)],t,squeeze(D(:,j,:)),squeeze(S(:,j,:)),objfn_constants.fn,objfn_constants.doses(j),objfn_constants.fn_opts),1:3)*objfn_constants.weights; % leave this here for now to remember the form of this function for the OxStudyFull 
F = @(p) arrayfun(@(j) rawError(objfn_constants.p_setup_fn(p),t,...
    D(j),objfn_constants.fn,C{j},objfn_constants.fn_opts),1:m)*objfn_constants.weights;

%% make sure pbest is best
[pbest,min_val] = fmincon(F,pbest,[],[],[],[],profile_params.lb,profile_params.ub,[],profile_params.opts);
pbest = reshape(pbest,[],1); % make sure it is a column vector
val_at_center = min_val;
npars = numel(pbest);

out = cell(npars,1);
for i = 1:npars
    if save_all_pars
        out{i} = [pbest;val_at_center];
    else
        out{i} = [pbest(i);val_at_center];
    end
    lb_temp = profile_params.lb;
    ub_temp = profile_params.ub;
    for dir = [-1,1] % move left and right along this parameter dimension
        if pbest(i)~=0
            dxi = pbest(i)*profile_params.initial_step_prop(i); % move at 1% of the best fit val
        else
            dxi = profile_params.smallest_par_step(i); % do not let the step size go below this as it steps towards the lower boundary
        end
        x0 = pbest;
        temp_val = zeros(1,profile_params.min_num_steps(i));
        if save_all_pars
            temp_par = zeros(npars,profile_params.min_num_steps(i));
            par_ind = i;
        else
            temp_par = zeros(1,profile_params.min_num_steps(i));
            par_ind = 1;
        end
        for j = 1:profile_params.min_num_steps(i)
            if x0(i) + dir*dxi > profile_params.para_ranges(i,2)
                dxi = dxi*profile_params.shrinking_factor^(ceil(log((profile_params.para_ranges(i,2)-x0(i))/dxi)/log(profile_params.shrinking_factor))); % if the step takes us below the lower threshold (usuallly 0), then take a smaller step (for the larger threshold, I have better reason to not let that grow too much)
            end
            x0(i) = x0(i) + dir*dxi;
            lb_temp(i) = x0(i);
            ub_temp(i) = x0(i);
            [x0,temp_val(j)] = fmincon(F,x0,[],[],[],[],lb_temp,ub_temp,[],profile_params.opts);
            min_val = min(min_val,temp_val(j));
            if save_all_pars
                temp_par(:,j) = x0;
            else
                temp_par(j) = x0(i);
            end
            if temp_val(j) >= min_val + profile_params.threshold
                x0(i) = x0(i) - dir*dxi; % reset to previous value to then increase by smaller step size
                dxi = dxi / pi(); % reduce step size, choose pi so we won't repeat any previously chosen values
            end
        end
        dxi = max(dxi,profile_params.min_par_step(i)); % for the remainder of the search, use this min value if dxi is too small
        if dir==-1
            max_steps = ceil((x0(i)-profile_params.para_ranges(i,1))/dxi);
        else
            max_steps = ceil((profile_params.para_ranges(i,2)-x0(i))/dxi);
        end
        extra_pars = zeros(size(temp_par,1),max_steps);
        extra_vals = zeros(1,max_steps);
        last_val = temp_val(end);
        j = 0;
        while (x0(i) + dir*dxi <= profile_params.para_ranges(i,2)) && (last_val < min_val + profile_params.threshold)
            if x0(i) + dir*dxi < profile_params.para_ranges(i,1)
                dxi = dxi*profile_params.shrinking_factor^(ceil(log((x0(i)-profile_params.para_ranges(i,1))/dxi)/log(profile_params.shrinking_factor))); % if the step takes us below the lower threshold (usuallly 0), then take a smaller step (for the larger threshold, I have better reason to not let that grow too much)
                if dxi < profile_params.smallest_par_step(i)
                    break;
                end
            end
            j = j+1;
            x0(i) = x0(i) + dir*dxi;
            if save_all_pars
                extra_pars(:,j) = x0;
            else
                extra_pars(j) = x0(i);
            end
            lb_temp(i) = x0(i);
            ub_temp(i) = x0(i);
            [x0,last_val] = fmincon(F,x0,[],[],[],[],lb_temp,ub_temp,[],profile_params.opts);
            extra_vals(j) = last_val;
            min_val = min(min_val,last_val);
        end
        out{i} = [out{i},[temp_par;temp_val],[extra_pars(:,1:j);extra_vals(1:j)]];
    end
    [~,order] = sort(out{i}(par_ind,:),"ascend");
    out{i} = out{i}(:,order);
end