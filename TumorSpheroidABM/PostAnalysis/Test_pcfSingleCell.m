clearvars;

% path_to_data_file = "../data/659808447356250/output_00000100.mat";
% 
% load(path_to_data_file,"tumor_locations","immune_locations","time")
% 
% tumor_locations = double(tumor_locations);
% immune_locations = double(immune_locations);
% 
% immune_ind = sub2ind([51,51,51],immune_locations(:,1),immune_locations(:,2),immune_locations(:,3));
% 
% all_subs = allCombos(1:51,1:51,1:51,'matlab');
% L = false([51,51,51]);
% L(immune_ind) = true;
% D = reshape(sqrt(sum(allCombos(0:50,0:50,0:50).^2,2)),[51,51,51]);
% 
% rr = linspace(0,26,91);
% 
% out = pcfSingleCell(tumor_locations(1,:),all_subs,L(:),D,rr,size(immune_locations,1)/numel(L),true);
% 
% %%
% yy = zeros(length(rr),1);
% figure;
% line_plot = plot(rr,yy);
% weights = zeros(length(rr),1);
% tumor_locations = tumor_locations(randperm(size(tumor_locations,1)),:);
% out = zeros(length(rr),500);
% out2 = zeros(length(rr),500);
% occupied = L(:);
% % for i = 1:500
% for i = 1:size(tumor_locations,1)
%     out(:,i) = pcfSingleCell(tumor_locations(i,:),all_subs,occupied,D,rr,size(immune_locations,1)/numel(L),true);
% 
% %     assert(all(out(:,i)==out2(:,i) | (isnan(out(:,i)) & isnan(out2(:,i)))));
% 
%     not_nan_log = ~isnan(out(:,i));
% %     line_plot.YData = (line_plot.YData * (i-1) + out(:,i)')/i;
%     line_plot.YData(not_nan_log) = (line_plot.YData(not_nan_log) .* weights(not_nan_log)' + out(not_nan_log,i)');
%     weights = weights + not_nan_log;
%     line_plot.YData(not_nan_log) = line_plot.YData(not_nan_log)./weights(not_nan_log)';
%     if mod(i,500)==0
%         title(sprintf("Using N = %d Tumor Cells",i))
%         drawnow
%     end
% end
% 
% %% 
% locs = tumor_locations(randperm(size(tumor_locations,1)),:);
% imm_locs = double(immune_locations);
% stop_at_edge = true;
% rr = linspace(0,26,27);
% out = pcfAllCells(locs,all_subs,imm_locs,rr,stop_at_edge);

%%

% rr = linspace(0,26,27);
% stop_at_edge = true;
% summarize = true;
% proportion_to_use = 0.01;
% min_to_use = 4e2;
% options = "high_antigen_to_immune";
% 
% out = pcfTimeSeries("../data/741615711346208",rr,stop_at_edge,summarize,proportion_to_use,min_to_use,options);

%%
addpath("~/Lattice-PCF-4-MATLAB/src")


%% prepare inputs
cohort_name = "cohort_221012222824126";
options = defaultOptions();
options.stop_at_edge = true;
options.is_cross_pcf = true;
options.multiple_points_per_site = false; % whether multiple points can occupy the same lattice site
options.proportion_to_use = 0.01;
options.min_to_use = 4E2;

options.rr = linspace(0,26,27);

summarize = true;
load_options = "low_antigen_to_immune";
out = pcfCohort(sprintf("../data/%s/%s.mat",cohort_name,cohort_name),options,summarize,load_options);

load(sprintf("../data/%s/%s.mat",cohort_name,cohort_name),"tracked")

out = reshape(out,size(tracked));
clear tracked

switch load_options
    case "default"
        save(sprintf("../data/%s/pcf",cohort_name),"-v7.3")
    case "tumor_to_active"
        save(sprintf("../data/%s/pcf_tumor_to_active",cohort_name),"-v7.3")
    case "high_antigen_to_immune"
        save(sprintf("../data/%s/pcf_high_antigen_to_immune",cohort_name),"-v7.3")
    case "low_antigen_to_immune"
        save(sprintf("../data/%s/pcf_low_antigen_to_immune",cohort_name),"-v7.3")

end