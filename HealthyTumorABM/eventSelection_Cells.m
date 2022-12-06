function [M,out] = eventSelection_Cells(M)

rate_matrix = computeRateMatrix(M);

% weighted count of cells in G0
M.tracked.healthy_types(M.i,1:2) = M.tracked.healthy_types(M.i,1:2) + sum(min(M.dt,M.healthy(:,M.I.proliferation_timer)))/M.dt * [1,-1];
M.tracked.tumor_types(M.i,1:2) = M.tracked.tumor_types(M.i,1:2) + sum(min(M.dt,M.tumor(:,M.I.proliferation_timer)))/M.dt * [1,-1];

prob_matrix.healthy = 1-exp(-rate_matrix.healthy.*[max(0,M.dt-M.healthy(:,M.I.proliferation_timer)),M.dt*ones(M.healthy_count,1)]); % use M.dt as time step except for proliferation, use M.dt - remaining time to wait before next proliferation
prob_matrix.tumor   = 1-exp(-rate_matrix.tumor  .*[max(0,M.dt-M.tumor(:,M.I.proliferation_timer))  ,M.dt*ones(M.tumor_count,1)]); % use M.dt as time step except for proliferation, use M.dt - remaining time to wait before next proliferation

event_log_array_healthy = rand(size(prob_matrix.healthy))<prob_matrix.healthy; % choose whether events happen
event_log_array_tumor = rand(size(prob_matrix.tumor))<prob_matrix.tumor; % choose whether events happen

inactive_log_healthy = all(~event_log_array_healthy,2); % index of all tumor cells that do nothing
[out.healthy.active_ind,out.healthy.events] = find(event_log_array_healthy); % index of all tumor cells that do something and the out.healthy.events they will do

inactive_log_tumor = all(~event_log_array_tumor,2); % index of all tumor cells that do nothing
[out.tumor.active_ind,out.tumor.events] = find(event_log_array_tumor); % index of all tumor cells that do something and the out.healthy.events they will do

final_event_healthy = zeros(M.healthy_count,1,'uint8');
final_event_healthy(inactive_log_healthy) = size(rate_matrix.healthy,2)+1; % all inactive tumor cells do nothing (which is marked by final event index + 1)

final_event_tumor = zeros(M.tumor_count,1,'uint8');
final_event_tumor(inactive_log_tumor) = size(rate_matrix.tumor,2)+1; % all inactive tumor cells do nothing (which is marked by final event index + 1)

u = rand(length(out.healthy.events),1); % random number to decide when event occurs
min_wait = (out.healthy.events==1).*M.healthy(out.healthy.active_ind,M.I.proliferation_timer); % minimum time waited until proliferation event can possibly occur
rates_healthy = rate_matrix.healthy(sub2ind(size(rate_matrix.healthy),out.healthy.active_ind,out.healthy.events)); % rates at which all the selected out.healthy.events happen
dts_healthy = M.dt - min_wait; % adjust total dts for proliferation

out.healthy.time_to_event = min_wait - log(1-u.*(1-exp(-rates_healthy.*dts_healthy)))./rates_healthy; % given that the event happens in this M.dt update, this determines when it occurs in that update step; this is why it's not just min_wait-log(rand())./rates

u = rand(length(out.tumor.events),1); % random number to decide when event occurs
min_wait = (out.tumor.events==1).*M.tumor(out.tumor.active_ind,M.I.proliferation_timer); % minimum time waited until proliferation event can possibly occur
rates_tumor = rate_matrix.tumor(sub2ind(size(rate_matrix.tumor),out.tumor.active_ind,out.tumor.events)); % rates at which all the selected out.tumor.events happen
dts_tumor = M.dt - min_wait; % adjust total dts for proliferation

out.tumor.time_to_event = min_wait - log(1-u.*(1-exp(-rates_tumor.*dts_tumor)))./rates_tumor; % given that the event happens in this M.dt update, this determines when it occurs in that update step; this is why it's not just min_wait-log(rand())./rates

[out.healthy.time_to_event,event_order] = sort(out.healthy.time_to_event);
out.healthy.active_ind = out.healthy.active_ind(event_order);
out.healthy.events = out.healthy.events(event_order);

[out.tumor.time_to_event,event_order] = sort(out.tumor.time_to_event);
out.tumor.active_ind = out.tumor.active_ind(event_order);
out.tumor.events = out.tumor.events(event_order);

%% remove proliferations that happen after apoptosis

apop_event_ind = find(out.healthy.events==2);

remove_ind = false(length(out.healthy.active_ind),1);
for i = 1:length(apop_event_ind)
    ai = apop_event_ind(i); % find event index corresponding to this apoptosis   
    remove_ind(ai+find(out.healthy.active_ind(ai+1:end)==out.healthy.active_ind(ai))) = true; % get ready to remove all out.healthy.events by this apoptotic cell that occur after apoptosis
end

out.healthy.active_ind(remove_ind) = []; % remove these post-apoptotic events
out.healthy.events(remove_ind) = []; % remove these post-apoptotic events
out.healthy.time_to_event(remove_ind) = []; % remove these post-apoptotic events

final_event_healthy(out.healthy.active_ind) = out.healthy.events; % this should set the event status to the latest-occuring event (really, this should mean that if both prolif and apop happen in that order, then the event will be labeled apop)
    
M.healthy(:,M.I.event) = final_event_healthy; % label the tumor based on the final event it performs

%% remove proliferations that happen after apoptosis

apop_event_ind = find(out.tumor.events==2);

remove_ind = false(length(out.tumor.active_ind),1);
for i = 1:length(apop_event_ind)
    ai = apop_event_ind(i); % find event index corresponding to this apoptosis   
    remove_ind(ai+find(out.tumor.active_ind(ai+1:end)==out.tumor.active_ind(ai))) = true; % get ready to remove all out.tumor.events by this apoptotic cell that occur after apoptosis
end

out.tumor.active_ind(remove_ind) = []; % remove these post-apoptotic events
out.tumor.events(remove_ind) = []; % remove these post-apoptotic events
out.tumor.time_to_event(remove_ind) = []; % remove these post-apoptotic events

final_event_tumor(out.tumor.active_ind) = out.tumor.events; % this should set the event status to the latest-occuring event (really, this should mean that if both prolif and apop happen in that order, then the event will be labeled apop)
    
M.tumor(:,M.I.event) = final_event_tumor; % label the tumor based on the final event it performs

%% interlace the two types
out.active_ind = [out.healthy.active_ind;out.tumor.active_ind];
out.events = [out.healthy.events;out.tumor.events];
out.time_to_event = [out.healthy.time_to_event;out.tumor.time_to_event];
out.cell_type = [repmat("healthy",length(out.healthy.active_ind),1);repmat("tumor",length(out.tumor.active_ind),1)];

[out.time_to_event,order] = sort(out.time_to_event);
out.active_ind = out.active_ind(order);
out.events = out.events(order);
out.cell_type = out.cell_type(order);


%% if plotting with prob_matrix.tumor, put that here
if M.plot_pars.plot_fig && mod(M.i,M.plot_pars.plot_every)==M.plot_pars.plot_offset
    plotFunction_EndEventSelection(M, prob_matrix);
end
