function [profile,pbest,val_at_center] = profileInOneDirection(profile,F,pbest,val_at_center,i,dir,profile_params,save_all_pars)

% pbest is the best fit at the point. the
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

npars = numel(pbest);
if dir==-1
    par_extremum = profile_params.para_ranges(i,1);
    par_exceeds_extremum = @(p) p < par_extremum;
    update_dxi = @(p,dxi) dxi*profile_params.shrinking_factor^(ceil(log((p-profile_params.para_ranges(i,1))/dxi)/log(profile_params.shrinking_factor))); % if the step takes us below the lower threshold (usuallly 0), then take a smaller step (for the larger threshold, I have better reason to not let that grow too much)
else
    par_extremum = profile_params.para_ranges(i,2);
    par_exceeds_extremum = @(p) p > par_extremum;
    update_dxi = @(p,dxi) dxi*profile_params.shrinking_factor^(ceil(log((profile_params.para_ranges(i,2)-p)/dxi)/log(profile_params.shrinking_factor))); % if the step takes us below the lower threshold (usuallly 0), then take a smaller step (for the larger threshold, I have better reason to not let that grow too much)
end

lb_temp = profile_params.lb;
ub_temp = profile_params.ub;
min_val = val_at_center;

if pbest(i)~=0
    dxi = abs(pbest(i))*profile_params.initial_step_prop(i); % start by moving at a proportion of the best value
else
    dxi = profile_params.smallest_par_step(i); % do not let the step size go below this as it steps towards a boundary
end
x0 = pbest;
first_vals = zeros(1,profile_params.min_num_steps(i));
if save_all_pars
    first_pars = zeros(npars,profile_params.min_num_steps(i));
else
    first_pars = zeros(1,profile_params.min_num_steps(i));
end
for j = 1:profile_params.min_num_steps(i)
    if (dir==-1 && x0(i)==profile_params.para_ranges(i,1)) || (dir==1 && x0(i)==profile_params.para_ranges(i,2))
        % then we're at the end of the range; so stop now
        first_pars = first_pars(:,1:j-1);
        first_vals = first_vals(1:j-1);
        if j==1
            return; % we started at the edge of the range so just be done and don't record anything
        else
            break;
        end
    end
    if par_exceeds_extremum(x0(i) + dir*dxi)
        dxi = update_dxi(x0(i),dxi);
    end
    x0_new = x0;
    x0_new(i) = x0_new(i) + dir*dxi;
    lb_temp(i) = x0_new(i);
    ub_temp(i) = x0_new(i);
    [x0_new,first_vals(j)] = fmincon(F,x0_new,profile_params.A,profile_params.b,[],[],lb_temp,ub_temp,[],profile_params.opts);
    if first_vals(j) < min_val
        min_val = first_vals(j);
        pbest = x0_new;
        val_at_center = min_val;
    end
    if save_all_pars
        first_pars(:,j) = x0_new;
    else
        first_pars(j) = x0_new(i);
    end
    if first_vals(j) >= min_val + profile_params.threshold
        % x0 = x0; reset to previous value to then increase by smaller step size
        dxi = dxi / pi(); % reduce step size, choose pi so we won't repeat any previously chosen values
    else
        x0 = x0_new; % continue on from the new x0 value
    end
end
dxi = dxi * profile_params.secondary_step_factor(i); % for the remainder of the search, use this value of dxi
if dir==-1
    max_steps = min(1000,ceil((x0(i)-profile_params.para_ranges(i,1))/dxi));
else
    max_steps = min(1000,ceil((profile_params.para_ranges(i,2)-x0(i))/dxi));
end
extra_pars = zeros(size(first_pars,1),max_steps);
extra_vals = zeros(1,max_steps);
last_val = first_vals(end);
j = 0;
max_allowable_step = Inf;
while true
    if x0(i)==par_extremum
        break; % already reached end of profile previously
    end
    if par_exceeds_extremum(x0(i) + dir*dxi)
        dxi = update_dxi(x0(i),dxi);
        if dxi < profile_params.smallest_par_step(i)
            dxi = abs(par_extremum - x0(i)); % force it to step exactly onto boundary
        end
    end
    x0_new = x0;
    x0_new(i) = x0_new(i) + dir*dxi;
    lb_temp(i) = x0_new(i);
    ub_temp(i) = x0_new(i);
    [x0_new,temp_val] = fmincon(F,x0_new,profile_params.A,profile_params.b,[],[],lb_temp,ub_temp,[],profile_params.opts);
    if log10(temp_val/last_val) > 1
        temp_val_prev = Inf;
        while log10(temp_val_prev/temp_val) > 1
            temp_val_prev = temp_val;
            [x0_new,temp_val] = fmincon(F,x0_new,profile_params.A,profile_params.b,[],[],lb_temp,ub_temp,[],profile_params.opts);
        end
    end
    if temp_val < min_val + profile_params.threshold
        j = j+1;
        if save_all_pars
            extra_pars(:,j) = x0_new;
        else
            extra_pars(j) = x0_new(i);
        end
        extra_vals(j) = temp_val;
        if temp_val < min_val
            min_val = temp_val;
            pbest = x0_new;
            val_at_center = min_val;
        end
        if dxi * profile_params.step_growth_factor(i) < max_allowable_step
            dxi = dxi * profile_params.step_growth_factor(i); % try taking bigger steps if staying within threshold
        else
            dxi = 0.5*(dxi+max_allowable_step);
        end
        x0 = x0_new;
        last_val = temp_val;
    elseif dxi <= profile_params.smallest_par_step(i) % then we'll assume that we're done; record last value and leave
        j = j+1;
        if save_all_pars
            extra_pars(:,j) = x0_new;
        else
            extra_pars(j) = x0_new(i);
        end
        extra_vals(j) = temp_val;
        break;
    else
        % return to previous step, shrink step size, and proceed
        max_allowable_step = dxi;
        if profile_params.step_growth_factor(i) > 1
            dxi = dxi * min(1 / profile_params.step_growth_factor(i), (min_val+profile_params.threshold-last_val)/(temp_val-last_val)); % linear interpolate to get step size that would land on threshold, but at least decrease by growth factor
        else
            dxi = 0.5 * dxi;
        end
    end

end
profile = [profile,[first_pars;first_vals],[extra_pars(:,1:j);extra_vals(1:j)]];