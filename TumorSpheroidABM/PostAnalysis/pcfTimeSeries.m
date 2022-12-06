function [out,t] = pcfTimeSeries(path_to_sim_folder,options,summarize,load_options)

addpath("~/Documents/MATLAB/myfunctions/")

warning("off",'MATLAB:load:variableNotFound')
load(sprintf("%s/output_constants.mat",path_to_sim_folder),"grid_size")
warning("on",'MATLAB:load:variableNotFound')

if ~exist("grid_size","var") % then the data was created before I started saving grid size
    load(sprintf("%s/output_constants.mat",path_to_sim_folder),"fgfr3_regions")
    grid_size = size(fgfr3_regions);
end

constants = prepareConstants(grid_size,options);

files = dir(sprintf("%s/output_*.mat",path_to_sim_folder));
files = files(arrayfun(@(x) matches(x.name,"output_" + digitsPattern(8) + (".mat")),files)); % only get numbered output files for analysis

t = zeros(1,numel(files));
if summarize
    out.avg = zeros(constants.nr,numel(files));
    out.std = zeros(constants.nr,numel(files));
else
    out = repmat(struct("G",[],"W",[]),[numel(files),1]);
end

for i = 1:numel(files)
    [t(i),subs.centers,subs.targets] = pcfLoadData([files(i).folder,'/',files(i).name],load_options);

    [G,W] = pcfAllCenters(subs,constants,options);
    if summarize
        out.avg(:,i) = sum(G.*W,2)./sum(W,2);
        for j = 1:constants.nr
            out.std(j,i) = std(G(j,:),W(j,:));
        end
    else
        out(i).G = G;
        out(i).W = W;
    end

%     if ~isempty(subs.targets)
%         pcfTimeSeriesPlot(constants.rr,t(1:i),out.avg(:,1:i)')
%         drawnow
%     end
end
        

end



% all_subs = allCombos(1:grid_size(1),1:grid_size(2),1:grid_size(3),'matlab');
% D = reshape(sqrt(sum(allCombos(0:grid_size(1)-1,0:grid_size(2)-1,0:grid_size(3)-1,'matlab').^2,2)),grid_size);
% 
% if isempty(fieldnames(buckets))
%     buckets.grid_size = grid_size;
%     if stop_at_edge
% 
%         center_subs = round(0.5*(grid_size+1));       
%         r_upper_bound = min([center_subs,grid_size-center_subs+1]); % minimum lattice steps to get one site off the lattice
%         rr = [rr(rr<r_upper_bound),r_upper_bound*(1-eps)];
%         d = sum((all_subs-center_subs).^2,2);
% 
%         buckets.volumes =  histcounts(d,rr.^2); % bring the outermost shell in a little so it does not get the closest "site" off the mesh
% 
% 
%     else
%         buckets.volumes = cell(grid_size);
%     end
%     bucket_ind = 1;
% else
%     bucket_idx_found = false;
%     for bucket_ind = 1:numel(buckets)
%         if isequal(buckets(bucket_ind).grid_size,grid_size) && ( (stop_at_edge && ~iscell(buckets(bucket_ind).volumes)) || (~stop_at_edge && iscell(buckets(bucket_ind).volumes)) )
%             bucket_idx_found = true;
%             break;
%         end
%     end
%     if ~bucket_idx_found
%         buckets(end+1).grid_size = grid_size;
%         if stop_at_edge
% 
%             center_subs = round(0.5*(grid_size+1));
%             d = sum((constants.all_subs-center_subs).^2,2);
% 
%             buckets.volumes =  histcounts(d,constants.rr.^2); % bring the outermost shell in a little so it does not get the closest "site" off the mesh
% 
% 
%         else
%             buckets.volumes = cell(grid_size);
%         end
%         bucket_ind = numel(buckets);
%     end
% end
% 
% 
% while exist(sprintf("%s/output_%08d.mat",path_to_sim_folder,fi),"file")~=0
%     [t(fi+1),source_locations,target_locations] = pcfLoadData(sprintf("%s/output_%08d.mat",path_to_sim_folder,fi),load_options);
%     if ~summarize
%         [out{fi+1},buckets(bucket_ind).volumes] = pcfAllCells(source_locations,target_locations,buckets(bucket_ind).volumes,grid_size,all_subs,D,rr,stop_at_edge,proportion_to_use,min_to_use,load_options);
%     else
%         [out_temp,buckets(bucket_ind).volumes] = pcfAllCells(source_locations,target_locations,buckets(bucket_ind).volumes,grid_size,all_subs,D,rr,stop_at_edge,proportion_to_use,min_to_use,load_options);
% 
%         out.avg(:,fi+1) = mean(out_temp,2,"omitnan");
%         out.sd(:,fi+1) = std(out_temp,[],2,"omitnan");
% 
%         pcfTimeSeriesPlot(rr,t,out.avg')
%     end
%     fi = fi+1;
% end