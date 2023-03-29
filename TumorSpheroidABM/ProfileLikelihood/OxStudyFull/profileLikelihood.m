function out = profileLikelihood(pbest,tt,data,data_std,C,phase_dependent_death,para_ranges,lb,ub,opts,threshold,min_par_step,save_all_pars)

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

F = @(p) sum(arrayfun(@(i) sum(((computeTimeSeries(p,tt,C(i),phase_dependent_death)-data(:,:,i))./data_std(:,:,i)).^2,'all'),1:numel(C)),'all');

%% make sure pbest is best
pbest = fmincon(F,pbest,[],[],[],[],lb,ub,[],opts);
pbest = reshape(pbest,[],1); % make sure it is a column vector
min_val = F(pbest);
if nargin < 13
    save_all_pars = false;
end
npars = numel(pbest);

out = cell(npars,1);
for i = 1:npars
    lb_temp = lb;
    ub_temp = ub;
    for dir = [-1,1] % move left and right along this parameter dimension
        dxi = pbest(i)*0.01; % move at 1% of the best fit val
        x0 = pbest;
        temp_val = zeros(1,10);
        if save_all_pars
            temp_par = zeros(npars,10);
        else
            temp_par = zeros(1,10);
        end
        for j = 1:10
            x0(i) = x0(i) + dir*dxi;
            lb_temp(i) = x0(i);
            ub_temp(i) = x0(i);
            [x0,temp_val(j)] = fmincon(F,x0,[],[],[],[],lb_temp,ub_temp,[],opts);
            if save_all_pars
                temp_par(:,j) = x0;
            else
                temp_par(j) = x0(i);
            end
            if temp_val(j) >= min_val + threshold
                x0(i) = x0(i) - dir*dxi; % reset to previous value to then increase by smaller step size
                dxi = dxi / pi(); % reduce step size, choose pi so we won't repeat any previously chosen values
            end
        end
        dxi = max(dxi,min_par_step(i)); % for the remainder of the search, use this min value
        while (x0(i) + dir*dxi >= para_ranges(i,1)) && (x0(i) + dir*dxi <= para_ranges(i,2)) && (temp_val(end) < min_val + threshold)
            x0(i) = x0(i) + dir*dxi;
            if save_all_pars
                temp_par(:,end+1) = x0;
            else
                temp_par(end+1) = x0(i);
            end
            lb_temp(i) = x0(i);
            ub_temp(i) = x0(i);
            [x0,temp_val(end+1)] = fmincon(F,x0,[],[],[],[],lb_temp,ub_temp,[],opts);
        end
        if dir==-1
            if save_all_pars
                [~,order] = sort(temp_par(i,:),'ascend'); % in case the algorithm went back and shrank dt
                out{i} = [temp_par(:,order),pbest;temp_val(order),min_val];
            else
                [temp_par,order] = sort(temp_par,'ascend');
                out{i} = [temp_par,pbest(i);temp_val(order),min_val];
            end
        else
            if save_all_pars
                [~,order] = sort(temp_par(i,:),'ascend'); % in case the algorithm went back and shrank dt
                out{i} = [out{i},[temp_par(order);temp_val(order)]];
            else
                [temp_par,order] = sort(temp_par,'ascend');
                out{i} = [out{i},[temp_par;temp_val(order)]];
            end
        end
    end

end