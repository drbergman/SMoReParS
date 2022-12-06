function M = initializeEvents(M)

M.events.times = [M.setup.set_tumor_time;M.setup.censor_date];
M.events.event_index = [1;Inf];

[M.events.times,order] = sort(M.events.times,1,"ascend");
M.events.event_index = M.events.event_index(order);

M.events.n = length(M.events.times);