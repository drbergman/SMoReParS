function out = profileLikelihood(pbest,objfn_data,profile_params,save_all_pars)

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

F = @(p) arrayfun(@(j) rawError(objfn_data.tt,objfn_data.count(:,j),objfn_data.count_std(:,j),objfn_data.state2_prop(:,j),objfn_data.state2_prop_std(:,j),p,objfn_data.fn,objfn_data.doses(j),objfn_data.fn_opts),1:3)*objfn_data.weights;

%% make sure pbest is best
pbest = fmincon(F,pbest,[],[],[],[],profile_params.lb,profile_params.ub,[],profile_params.opts);
pbest = reshape(pbest,[],1); % make sure it is a column vector
min_val = F(pbest);
npars = numel(pbest);

out = cell(npars,1);
for i = 1:npars
    lb_temp = profile_params.lb;
    ub_temp = profile_params.ub;
    for dir = [-1,1] % move left and right along this parameter dimension
        if pbest(i)~=0
            dxi = pbest(i)*profile_params.initial_step_prop(i); % move at 1% of the best fit val
        else
            dxi = 1e-4;
        end
        x0 = pbest;
        temp_val = zeros(1,profile_params.min_num_steps(i));
        if save_all_pars
            temp_par = zeros(npars,profile_params.min_num_steps(i));
        else
            temp_par = zeros(1,profile_params.min_num_steps(i));
        end
        for j = 1:profile_params.min_num_steps(i)
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
        while (x0(i) + dir*dxi >= profile_params.para_ranges(i,1)) && (x0(i) + dir*dxi <= profile_params.para_ranges(i,2)) && (temp_val(end) < min_val + profile_params.threshold)
            x0(i) = x0(i) + dir*dxi;
            if save_all_pars
                temp_par(:,end+1) = x0;
            else
                temp_par(end+1) = x0(i);
            end
            lb_temp(i) = x0(i);
            ub_temp(i) = x0(i);
            [x0,temp_val(end+1)] = fmincon(F,x0,[],[],[],[],lb_temp,ub_temp,[],profile_params.opts);
            min_val = min(min_val,temp_val(end));
        end
        if dir==-1
            if save_all_pars
                [~,order] = sort(temp_par(i,:),'ascend'); % in case the algorithm went back and shrank dxi
                out{i} = [temp_par(:,order),pbest;temp_val(order),min_val];
            else
                [temp_par,order] = sort(temp_par,'ascend');
                out{i} = [temp_par,pbest(i);temp_val(order),min_val];
            end
        else
            if save_all_pars
                [~,order] = sort(temp_par(i,:),'ascend'); % in case the algorithm went back and shrank dxi
                out{i} = [out{i},[temp_par(order);temp_val(order)]];
            else
                [temp_par,order] = sort(temp_par,'ascend');
                out{i} = [out{i},[temp_par;temp_val(order)]];
            end
        end
    end

end