function simCohort(M,cohort_pars)

nsamps_per_condition = cohort_pars.nsamps_per_condition;

% lattice sampling of fields of M
fn = string(fieldnames(M));
cohort.lattice_parameters = struct("path",{},"values",{});
for i = 1:numel(fn)
    current_struct_path = fn(i);
    cohort.lattice_parameters = grabFields(M.(fn(i)),cohort.lattice_parameters,current_struct_path);
end

if cohort_pars.link_arrest_coeffs % enforce that the arrest coefficients are identical
    for li = 1:size(cohort_pars.linkings,1) % linking index
        arrest_coeff_inds = [];
        vals = {};
        paths = {};
        for i = 1:numel(cohort.lattice_parameters)
            if ~iscell(cohort.lattice_parameters(i).path) && startsWith(cohort.lattice_parameters(i).path(end),"arrest_coeff") && any(endsWith(cohort.lattice_parameters(i).path(end),cohort_pars.linkings(li,:)))
                arrest_coeff_inds(end+1) = i;
                vals{end+1} = cohort.lattice_parameters(i).values;
                paths{end+1} = cohort.lattice_parameters(i).path;
            end
        end
        if length(arrest_coeff_inds)>1 % only make this change if two arrest coeffs are varied
            assert(numel(unique(cellfun(@numel,vals)))==1) % make sure that each of these values has the same number of elements
            cohort.lattice_parameters(end+1).values = cat(2,vals{:});
            cohort.lattice_parameters(end).path = paths;
            cohort.lattice_parameters(arrest_coeff_inds) = [];
        end
    end
end

cohort.all_parameters = flattenStruct(M);
all_fn = fieldnames(cohort.all_parameters);


cohort_size = arrayfun(@(i) size(cohort.lattice_parameters(i).values,1),1:numel(cohort.lattice_parameters));
total_runs = prod(cohort_size) * nsamps_per_condition;
cohort_pars.total_runs = total_runs;

colons = repmat({':'},[1,length(cohort_size)]);
vp_ind = cell(1,length(cohort_size));

% sims_to_check = dir("data/sims/*");
% sims_to_check = sims_to_check([sims_to_check.isdir]);
% for i = numel(sims_to_check):-1:1
%     if ~startsWith(sims_to_check(i).name,digitsPattern(1)) % I only want folders that start with a number (not "." or ".." etc)
%         sims_to_check(i) = [];
%     else
%         sims_to_check(i).name = string(sims_to_check(i).name);
%     end
% end
% sims_to_check = [sims_to_check.name];
cohort.ids = repmat("",[cohort_size,nsamps_per_condition,1]); % put the 1 at the end in case cohort_size = []; this way it creates a 1D vector of ids rather than a square array of ids...silly matlab
n_found = 0;


%% check if previous cohorts ran these sims
previous_cohorts = dir("data/cohort_*");
for i = 1:numel(previous_cohorts)

    PC = load(sprintf("data/%s/output.mat",previous_cohorts(i).name));
    % sims_to_check = setdiff(sims_to_check,PC.ids(:));
    if ~isequal(sort(string(all_fn)),sort(string(fieldnames(PC.all_parameters))))
        % these did not have the same parameters coming in, so move on
        continue;
    end

    skip_this_cohort = false;
    D = cell(length(cohort_size),length(PC.lattice_parameters)); % tracks which indices in each dim pair should be mapped
    current_fixed = cell(length(cohort_size),1); % the possible indices in each dimension that have been found in current
    previous_fixed = cell(length(PC.lattice_parameters),1); % the possible indices in each dimension that have been found in previous
    current_fixed_inds = cell(length(cohort_size),1); % the indices in each dimension of the cohort that have some samples found, but the previous cohort did not vary along this dimension
    previous_inds = cell(length(PC.lattice_parameters),1); % the indices in each dimension of the PC that will be mapped to the new ids
    dim_match = false(length(cohort_size),length(PC.lattice_parameters)); % a matrix accounting for which current dims get mapped to from previous dims
    %     dim_match = zeros(0,2); % match target dims with copy dims (in case they are in different orders)

    for j = 1:numel(all_fn)

        current_val = cohort.all_parameters.(all_fn{j});
        previous_val = PC.all_parameters.(all_fn{j});
        if size(current_val,1)==1 % then this parameter is not currently being varied, check that at least one of previous used it
            if size(previous_val,1)==1 % then this was also not varied previously
                if ~isequal(current_val,previous_val)
                    skip_this_cohort = true;
                    break;
                end
            else % then need this was varied previously, need to see if the current val matches any of the previous, and then record where those are
                prev_dim = NaN; % clear any value from here
                [used_prev,prev_ind] = ismember(current_val,previous_val,"rows");
                if ~used_prev % make sure the current value was one of the previously used
                    skip_this_cohort = true;
                    break;
                end
                for k = 1:numel(PC.lattice_parameters)
                    if iscell(PC.lattice_parameters(k).path) % then this lattice parameter was a combination of parameters
                        lattice_par_found = false;
                        for l = 1:numel(PC.lattice_parameters(k).path)
                            if strcmp(all_fn{j},[PC.lattice_parameters(k).path{l}{1},'_DOT_',PC.lattice_parameters(k).path{l}{2}])
                                prev_dim = k;
                                lattice_par_found = true;
                                break;
                            end
                        end
                        if lattice_par_found
                            break; % stop looping through PC lattice parameters
                        end
                    else
                        if strcmp(all_fn{j},[PC.lattice_parameters(k).path{1},'_DOT_',PC.lattice_parameters(k).path{2}])
                            prev_dim = k;
                            break;
                        end
                    end
                end

                [previous_fixed,skip_this_cohort] = updateFixed(previous_fixed,prev_dim,prev_ind,skip_this_cohort);
                if skip_this_cohort
                    break;
                end


            end

        else % then this parameter is currently being varied

            current_dim = NaN; % clear any value from here
            if size(previous_val,1)==1 % then this was not varied previously
                [re_use,current_ind] = ismember(previous_val,current_val,"rows");
                if ~re_use % then none of these were used before, move on
                    skip_this_cohort = true;
                    break;
                end
                for k = 1:numel(cohort.lattice_parameters)
                    if iscell(cohort.lattice_parameters(k).path) % then this lattice parameter was a combination of parameters
                        lattice_par_found = false;
                        for l = 1:numel(cohort.lattice_parameters(k).path)
                            if strcmp(all_fn{j},[cohort.lattice_parameters(k).path{l}{1},'_DOT_',cohort.lattice_parameters(k).path{l}{2}])
                                current_dim = k;
                                lattice_par_found = true;
                                break;
                            end
                        end
                        if lattice_par_found
                            break; % stop looping through cohort lattice parameters
                        end
                    else
                        if strcmp(all_fn{j},[cohort.lattice_parameters(k).path{1},'_DOT_',cohort.lattice_parameters(k).path{2}])
                            current_dim = k;
                            break;
                        end
                    end
                end

                [current_fixed,skip_this_cohort] = updateFixed(current_fixed,current_dim,current_ind,skip_this_cohort);
                if skip_this_cohort
                    break;
                end

            else % it was varied previously as well

                prev_dim = NaN; % clear any value from here
                for k = 1:numel(PC.lattice_parameters)
                    if iscell(PC.lattice_parameters(k).path) % then this lattice parameter was a combination of parameters
                        lattice_par_found = false;
                        for l = 1:numel(PC.lattice_parameters(k).path)
                            if strcmp(all_fn{j},[PC.lattice_parameters(k).path{l}{1},'_DOT_',PC.lattice_parameters(k).path{l}{2}])
                                prev_dim = k;
                                lattice_par_found = true;
                                break;
                            end
                        end
                        if lattice_par_found
                            break; % stop looping through PC lattice parameters
                        end
                    else
                        if strcmp(all_fn{j},[PC.lattice_parameters(k).path{1},'_DOT_',PC.lattice_parameters(k).path{2}])
                            prev_dim = k;
                            break;
                        end
                    end
                end

                for k = 1:numel(cohort.lattice_parameters)
                    if iscell(cohort.lattice_parameters(k).path) % then this lattice parameter was a combination of parameters
                        lattice_par_found = false;
                        for l = 1:numel(cohort.lattice_parameters(k).path)
                            if strcmp(all_fn{j},[cohort.lattice_parameters(k).path{l}{1},'_DOT_',cohort.lattice_parameters(k).path{l}{2}])
                                current_dim = k;
                                lattice_par_found = true;
                                break;
                            end
                        end
                        if lattice_par_found
                            break; % stop looping through cohort lattice parameters
                        end
                    else
                        if strcmp(all_fn{j},[cohort.lattice_parameters(k).path{1},'_DOT_',cohort.lattice_parameters(k).path{2}])
                            current_dim = k;
                            break;
                        end
                    end
                end

                match_found = false;
                current_already_fixed = ~isempty(current_fixed{current_dim});
                previous_already_fixed = ~isempty(previous_fixed{prev_dim});
                already_mapped = ~isempty(D{current_dim,prev_dim});
                if previous_already_fixed
                    prev_iter = previous_fixed{prev_dim};
                else
                    prev_iter = 1:size(previous_val,1);
                end
                for k = prev_iter
                    [re_use,current_ind] = ismember(previous_val(k,:),current_val,"rows");
                    if re_use && (~current_already_fixed || ismember(current_ind,current_fixed{current_dim})) && (~already_mapped || ismember([current_ind,k],D{current_dim,prev_dim},"rows"))
                        match_found = true;
                        if ~already_mapped
                            D{current_dim,prev_dim}(end+1,:) = [current_ind,k];
                        end
                        if ~current_already_fixed
                            current_fixed{current_dim}(end+1) = current_ind;
                        end
                        if ~previous_already_fixed
                            previous_fixed{prev_dim}(end+1) = k;
                        end
                    elseif previous_already_fixed % if this ind cannot be found, make sure to remove it from previous_fixed
                        previous_fixed{prev_dim} = setdiff(previous_fixed{prev_dim},k);
                        [D,current_fixed,previous_fixed] = scanForDeletions(D,current_fixed,previous_fixed,true,prev_dim,k);
                    end
                end
                if ~match_found
                    skip_this_cohort = true;
                    break;
                end
            end
        end
    end

    if skip_this_cohort
        continue;
    end

    target_sz = cellfun(@numel,current_fixed);

    target_inds = cell(length(target_sz),1);
    ti = cell(length(target_sz),1);
    for j = 1:prod(target_sz)
        [target_inds{:}] = ind2sub(target_sz,j);
        for k = 1:length(target_sz)
            ti{k} = current_fixed{k}(target_inds{k});
        end
        ci = matchToCopy(ti,D,previous_fixed);
        if isempty(ci)
            continue;
        end
        temp_ids = cohort.ids(ti{:},:);
        blank_log = temp_ids=="";
        n_blank = sum(blank_log);
        if n_blank==0 % then all have been found, move on
            continue;
        end
        new_ids_temp = PC.ids(ci{:},:);
        new_ids_temp = setdiff(new_ids_temp,temp_ids);
        blanks_to_fill = find(blank_log);
        if numel(new_ids_temp) > n_blank
            new_ids_temp = new_ids_temp(randperm(numel(new_ids_temp),n_blank)); % select random ids from before
        elseif numel(new_ids_temp) < n_blank
            blanks_to_fill = blanks_to_fill(1:numel(new_ids_temp));
        end
        temp_ids(blanks_to_fill) = new_ids_temp;
        cohort.ids(ti{:},:) = sort(temp_ids);
        n_found = n_found + numel(blanks_to_fill);
    end

    if n_found == total_runs %|| isempty(sims_to_check)
        break;
    end

end

if ~isempty(cohort_size) % if its empty, then cohort.ids should just be a column vector already, as desired
    cohort.ids = squeeze(permute(cohort.ids,[length(cohort_size)+1,1:length(cohort_size)])); % put the sample dimension along first dimension for sims
end

% sims_to_check = sims_to_check(randperm(numel(sims_to_check))); % to not bias samples towards the first sims I ran

% %% record which previous sims have the right parameters
% fn = fieldnames(M);
% for i = 1:numel(sims_to_check)
%     if exist(sprintf("data/sims/%s/output_constants.mat",sims_to_check(i)),"file") && exist(sprintf("data/sims/%s/output_final.mat",sims_to_check(i)),"file")
%         X = load(sprintf("data/sims/%s/output_constants.mat",sims_to_check(i)));
%         these_match = true;
%         for j = 1:numel(fn)
%             if ~strcmp(fn{j},"plot_pars") % don't worry about plot_pars being equal
%                 if ~isfield(X,fn{j}) % then this wasn't a saved field before, so don't consider this sim
%                     these_match = false;
%                     break;
%                 end
%                 par_fn = fieldnames(M.(fn{j}));
%                 for k = 1:numel(par_fn)
%                     if ~isfield(X.(fn{j}),par_fn{k}) % then this wasn't a saved field before, so don't consider this sim
%                         these_match = false;
%                         break;
%                     end
%                     if strcmp(par_fn{k},'n_regions') % n_regions is actually set for each substrate (i.e. I can remove n_regions from base parameter stuff)
%                         continue;
%                     end
%                     if strcmp(par_fn{k},'deactivation_function') % it seems the isequal cannot really check if two anonymous functions are equal
%                         continue;
%                     end
%                     if startsWith(par_fn{k},'grid_size_microns_') && M.setup.use_carrying_capacity_for_grid_size
%                         continue; % if using the carrying capacity, then the grid size does not matter at this point
%                     end
%                     if size(M.(fn{j}).(par_fn{k}),1)==1 % then this parameter is not being varied
%                         if ~isequal(M.(fn{j}).(par_fn{k}),X.(fn{j}).(par_fn{k}))
%                             these_match = false;
%                             break;
%                         end
%                     end
%                 end
%                 if ~these_match
%                     break;
%                 end
%             end
%         end
%         if ~these_match % then already know something doesn't match, don't need to compare values
%             break;
%         end
%         for vpi = 1:numel(cohort.lattice_parameters)
%             xtemp = X;
%             for pi = 1:length(cohort.lattice_parameters(vpi).path)
%                 xtemp = xtemp.(cohort.lattice_parameters(vpi).path{pi});
%             end
%             vp_ind{vpi} = find(cohort.lattice_parameters(vpi).values==xtemp);
%             if isempty(vp_ind{vpi})
%                 these_match = false;
%                 break;
%             end
%         end
%         if these_match
%             sample_ind = find(cohort.ids(:,vp_ind{:})=="",1);
%             if ~isempty(sample_ind) && sample_ind <= nsamps_per_condition
%                 cohort.ids(sample_ind,vp_ind{:}) = sims_to_check(i);
%                 n_found = n_found+1;
%                 if n_found>=total_runs
%                     break;
%                 end
%             end
% 
%         end
%     end
% end

inds_to_run = find(cohort.ids=="");
total_runs = numel(inds_to_run);

%% now fill out the rest of the sim array/grab the sim data identified above
cohort.ids = cohort.ids(:);
cohort_start_time = string(datetime("now","Format","yyMMddHHmm"));
if ~isfield(cohort_pars,"cohort_identifier")
    cohort_pars.cohort_identifier = cohort_start_time; % default to this for determining an id if none given
    while exist(sprintf("data/cohort_%s",cohort_pars.cohort_identifier),"dir") % just in case this directory already exists somehow (not sure how to processes could start at the same time to the millisecond and then one create this folder before the other looks for it)
        cohort_pars.cohort_identifier = string(datetime("now","Format","yyMMddHHmmss")); % default to this for determining an id if none given
    end
end


if total_runs>=cohort_pars.min_parfor_num
    F(1:total_runs) = parallel.FevalFuture;
    ppool = gcp;
    cohort_pars.num_workers = ppool.NumWorkers;
    for ri = 1:total_runs % run index
        if mod(ri,1000)==0
            fprintf("Setting up simulation %d of %d...\n",ri,total_runs)
        end
        [~,vp_ind{colons{:}}] = ind2sub([nsamps_per_condition,cohort_size],inds_to_run(ri));
        for vpi = 1:numel(cohort.lattice_parameters)
            M = setField(M,cohort.lattice_parameters(vpi).path,cohort.lattice_parameters(vpi).values(vp_ind{vpi},:));
        end
        M.save_pars.sim_identifier = sprintf("%s_%d",cohort_pars.cohort_identifier,ri);
        F(ri) = parfeval(ppool,@simPatient,1,M);
    end
else
    cohort_pars.num_workers = 1;
end

%%
cohort_pars.mu_n = 0;
cohort_pars.start = tic;
cohort_pars.batch_start = tic;
cohort_pars.total_runs = total_runs; % set this for calculating time remaining

for ri = total_runs:-1:1
    if total_runs>=cohort_pars.min_parfor_num
        [temp,out_temp] = fetchNext(F);
        idx = inds_to_run(temp);
    else
        idx = inds_to_run(ri);
        [~,vp_ind{colons{:}}] = ind2sub([nsamps_per_condition,cohort_size],idx);
        for vpi = 1:numel(cohort.lattice_parameters)
            M = setField(M,cohort.lattice_parameters(vpi).path,cohort.lattice_parameters(vpi).values(vp_ind{vpi},:));
        end
        M.save_pars.sim_identifier = sprintf("%s_%d",cohort_pars.cohort_identifier,ri);
        out_temp = simPatient(M);
    end
    cohort_pars = updateCohortTimer(cohort_pars,total_runs-ri+1);
    if isfield(out_temp.save_pars,"sim_identifier")
        cohort.ids(idx) = out_temp.save_pars.sim_identifier;
    end
end

if M.save_pars.make_save < Inf
    cohort.ids = reshape(cohort.ids,[nsamps_per_condition,cohort_size,1]);
end

if ~isempty(cohort_size)
    cohort.ids = squeeze(permute(cohort.ids,[2:length(cohort_size)+1,1])); % put the sample dimension back along last dimension
end

mkdir(sprintf("data/cohort_%s",cohort_pars.cohort_identifier))

save(sprintf("data/cohort_%s/output",cohort_pars.cohort_identifier),"nsamps_per_condition","total_runs","cohort_size")
save(sprintf("data/cohort_%s/output",cohort_pars.cohort_identifier),'-struct',"cohort","-append")

fprintf("Finished cohort. Folder is: cohort_%s\n",cohort_pars.cohort_identifier)

if cohort_pars.check_cohort_grab
    checkCohortGrab(cohort_pars.cohort_identifier)
end

end

function lattice_parameters = grabFields(S,lattice_parameters,incoming_struct_path)

fn = string(fieldnames(S));
for i = 1:numel(fn)
    current_struct_path = [incoming_struct_path,fn(i)];
    if isstruct(S.(fn(i)))
        lattice_parameters = grabFields(S.(fn(i)),lattice_parameters,current_struct_path);
    elseif size(S.(fn(i)),1)>1 % then vary over these parameters
        lattice_parameters(end+1) = struct("path",current_struct_path,"values",S.(fn(i)));
    end
end

end

function [sim_this,sims_to_check,new_start_ind,sim_id] = findSimilarSims(M,sims_to_check,start_ind)
sim_id = "";
new_start_ind = numel(sims_to_check)+1; % if no match is found for these pars, then don't search for a match next time
sim_this = true;
for i = start_ind:numel(sims_to_check) % look for the first one that matches the inputs here
    if exist(sprintf("data/sims/%s/output_constants.mat",sims_to_check(i)),"file") && exist(sprintf("data/sims/%s/output_final.mat",sims_to_check(i)),"file")
        X = load(sprintf("data/sims/%s/output_constants.mat",sims_to_check(i)));
        fn = fieldnames(M);
        these_match = true;
        for j = 1:numel(fn)
            if (~strcmp(fn{j},"plot_pars")) % don't worry about plot_pars being equal
                par_fn = fieldnames(M.(fn{j}));
                for k = 1:numel(par_fn)
                    if strcmp(par_fn{k},'n_regions') % n_regions is actually set for each substrate (i.e. I can remove n_regions from base parameter stuff)
                        continue;
                    end
                    if strcmp(par_fn{k},'deactivation_function') % it seems the isequal cannot really check if two anonymous functions are equal
                        continue;
                    end
                    if startsWith(par_fn{k},'grid_size_microns_') && M.setup.use_carrying_capacity_for_grid_size
                        continue; % if using the carrying capacity, then the grid size does not matter at this point
                    end
                    if ~isfield(X,fn{j}) || ~isfield(X.(fn{j}),par_fn{k}) || ~isequal(M.(fn{j}).(par_fn{k}),X.(fn{j}).(par_fn{k}))
                        these_match = false;
                        break;
                    end
                end
                if ~these_match
                    break;
                end
            end
        end
        if these_match
            sim_this = false;
            sim_id = string(sims_to_check(i));
            sims_to_check(i) = [];
            new_start_ind = i;
            return;
        end
    end
end
end

function out = grabSimData(folder_name)

out.save_pars.sim_identifier = folder_name;

end

function ids = reuseSims(ids,copy_ids,sample_dim)

if sample_dim == 1 % then these are both vectors of IDs
    n = length(ids);
    ids = [ids;copy_ids(copy_ids~="")];
    ids = ids(1:n);
else
    for i = 1:size(ids,1)
        colons = repmat(':',1,sample_dim-1);
        ids(i,colons{:}) = reuseSims(squeeze(ids(i,colons{:})),squeeze(copy_ids(i,colons{:})),sample_dim-1);
    end
end
end


function [fixed,skip_this_cohort] = updateFixed(fixed,dim,ind,skip_this_cohort)

if isempty(fixed{dim})
    fixed{dim} = ind;
else
    fixed{dim} = intersect(fixed{dim},ind);
    if isempty(fixed{dim})
        skip_this_cohort = true; % then this is in conflict with the previously found indices for this dimension
    end
end
end

function copy_inds = matchToCopy(ti,D,previous_fixed)

n = length(previous_fixed);
copy_inds = cell(n,1);
for k = 1:n
    if numel(previous_fixed{k})==1
        copy_inds{k} = previous_fixed{k}; % only one to grab here, so just grab it
    else
        ind = NaN;
        for ri = 1:size(D,1) % loop over target dims to find the ones that match with dim k in previous
            if ~isempty(D{ri,k})
                if isnan(ind)
                    ind = D{ri,k}((D{ri,k}(:,1)==ti{ri}),2);
                    assert(numel(ind)==1) % should only find a single match here
                elseif ~isequal(ind,D{ri,k}((D{ri,k}(:,1)==ti{ri}),2)) % should match the previous index found; skip this one
                    copy_inds = [];
                    return;
                end
            end
        end
        if isnan(ind) % could not find the proper previous indices to match the target indices; I don't think this should ever happen
            error("No index found???")
            %             copy_inds = [];
            %             return;
        end
        copy_inds{k} = ind;
    end
end
end

function [D_full,current_fixed,previous_fixed] = scanForDeletions(D_full,current_fixed,previous_fixed,fix_previous,dim,I)

if fix_previous
    D = D_full(:,dim);
    col_ind = 2;
    other_ind = 1;
else
    D = D_full(dim,:);
    col_ind = 1;
    other_ind = 2;
end

for i = 1:length(D) % find where it was previously matched and remove those
    if ~isempty(D{i})
        ind = D{i}(:,col_ind) == I;
        if ~any(ind)
            continue;
        end
        other_to_remove = D{i}(ind,other_ind);
        if fix_previous
            D_full{i,dim}(ind,:) = [];
            current_fixed{i} = setdiff(current_fixed{i},other_to_remove);
            [D_full,current_fixed,previous_fixed] = scanForDeletions(D_full,current_fixed,previous_fixed,false,i,other_to_remove);
        else
            D_full{dim,i}(ind,:) = [];
            previous_fixed{i} = setdiff(previous_fixed{i},other_to_remove);
            [D_full,current_fixed,previous_fixed] = scanForDeletions(D_full,current_fixed,previous_fixed,true,i,other_to_remove);
        end
    end
end


end
