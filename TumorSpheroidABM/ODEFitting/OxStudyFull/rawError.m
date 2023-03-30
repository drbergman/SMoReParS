function out = rawError(tt,y,s,prop,s_prop,p,fn,dose,fn_opts)

sim_data = fn(p,tt,dose,fn_opts);
% sim_total = sum(sim_data,2);
% 
% if size(sim_data,2)==4
%     sim_prop = sum(sim_data(:,3:4),2)./sim_total;
% else
%     sim_prop = sim_data(:,2)./sim_total;
% end
% 
% out = sum(((sim_total-y)./s).^2,"all","omitnan") + sum(((sim_prop-prop)./s_prop).^2,"all","omitnan");

out = sum(((sim_data - [y,prop])./[s,s_prop]).^2,"all","omitnan");

disp('')
