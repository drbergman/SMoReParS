function simCohort(M,cohort_pars)

nsamps_per_condition = cohort_pars.nsamps_per_condition;

% lattice sampling of fields of M
fn = string(fieldnames(M));
cohort.lattice_parameters = struct("path",{},"values",{});
for i = 1:numel(fn)
    current_struct_path = fn(i);
    cohort.lattice_parameters = grabFields(M.(fn(i)),cohort.lattice_parameters,current_struct_path);
end

% non-lattice parameter varying (for dosing schedule, rectangle sizes...and others?)
if cohort_pars.last_dose_is_no_dose
    dose_start_inds = [];
    vals = {};
    paths = {};
    for i = 1:numel(cohort.lattice_parameters)
        if cohort.lattice_parameters(i).path(end)=="start_day"
            dose_start_inds(end+1) = i;
            vals{end+1} = cohort.lattice_parameters(i).values;
            paths{end+1} = cohort.lattice_parameters(i).path;
        end
    end
    cohort.lattice_parameters(end+1).values = allCombos(vals{:},'matlab');
    cohort.lattice_parameters(end).values(end,:) = Inf;
    cohort.lattice_parameters(end).path = paths;
    cohort.lattice_parameters(dose_start_inds) = [];
end
cohort_size = arrayfun(@(i) size(cohort.lattice_parameters(i).values,1),1:numel(cohort.lattice_parameters));
total_runs = prod(cohort_size) * nsamps_per_condition;
cohort_pars.total_runs = total_runs;

colons = repmat({':'},[1,length(cohort_size)]);
vp_ind = cell(1,length(cohort_size));

sims_to_check = dir("data/*");
sims_to_check = sims_to_check([sims_to_check.isdir]);
for i = numel(sims_to_check):-1:1
    if ~startsWith(sims_to_check(i).name,digitsPattern(1)) % I only want folders that start with a number 
        sims_to_check(i) = [];
    end
end

sims_to_check = sims_to_check(randperm(numel(sims_to_check))); % to not bias samples towards the first sims I ran
cohort.ids = repmat("",[nsamps_per_condition,cohort_size]);

%% record which previous sims have the right parameters
n_found = 0;
fn = fieldnames(M);
for i = 1:numel(sims_to_check)
    if exist(sprintf("data/%s/output_constants.mat",sims_to_check(i).name),"file") && exist(sprintf("data/%s/output_final.mat",sims_to_check(i).name),"file")
        X = load(sprintf("data/%s/output_constants.mat",sims_to_check(i).name));
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
                    if size(M.(fn{j}).(par_fn{k}),1)==1 % then this parameter is being varied
                        if ~isfield(X,fn{j}) || ~isfield(X.(fn{j}),par_fn{k}) || ~isequal(M.(fn{j}).(par_fn{k}),X.(fn{j}).(par_fn{k}))
                            these_match = false;
                            break;
                        end
                    end
                end
                if ~these_match
                    break;
                end
            end
        end
        for vpi = 1:numel(cohort.lattice_parameters)
            xtemp = X;
            for pi = 1:length(cohort.lattice_parameters(vpi).path)
                xtemp = xtemp.(cohort.lattice_parameters(vpi).path{pi});
            end
            vp_ind{vpi} = find(cohort.lattice_parameters(vpi).values==xtemp);
            if isempty(vp_ind{vpi})
                these_match = false;
                break;
            end
        end
        if these_match
            sample_ind = find(cohort.ids(:,vp_ind{:})=="",1);
            if ~isempty(sample_ind) && sample_ind <= nsamps_per_condition
                cohort.ids(sample_ind,vp_ind{:}) = sims_to_check(i).name;
                n_found = n_found+1;
                if n_found>=total_runs
                    break;
                end
            end
            
        end
    end
end

sims_to_check(:) = [];

%% now fill out the rest of the sim array/grab the sim data identified above
cohort.ids = cohort.ids(:);

ids = reshape(cohort.ids,[nsamps_per_condition,cohort_size]); % use this to make sure we don't use a simulation 2x

if total_runs>=cohort_pars.min_parfor_num
    F(1:total_runs) = parallel.FevalFuture;
    ppool = gcp;
    cohort_pars.num_workers = ppool.NumWorkers;
    for ri = 1:total_runs % run index
        fprintf("Setting up simulation #%d...\n",ri)
        [sample_ind,vp_ind{colons{:}}] = ind2sub([nsamps_per_condition,cohort_size],ri);
        if sample_ind == 1
            start_ind = 1;
        end
        if cohort.ids(ri)~="" % if the above code found a sim here, then grab that
            F(ri) = parfeval(ppool,@grabSimData,1,cohort.ids(ri));
        else
            for vpi = 1:numel(cohort.lattice_parameters)
                M = setField(M,cohort.lattice_parameters(vpi).path,cohort.lattice_parameters(vpi).values(vp_ind{vpi},:));
            end
            [sim_this,sims_to_check,start_ind,sim_id] = findSimilarSims(M,sims_to_check,start_ind);
            if sim_this || any(ids(:,vp_ind{:})==sim_id) % make sure that we didn't already record this sim for these parameter values
                F(ri) = parfeval(ppool,@simPatient,1,M);
            else % if the above code did not (somehow) grab this sim, then grab this sim data now
                F(ri) = parfeval(ppool,@grabSimData,1,sim_id);
            end
        end
    end
else
    cohort_pars.num_workers = 1;
end

%% 
cohort_pars.mu_n = 0;
cohort_pars.start = tic;
cohort_pars.batch_start = tic;

for ri = total_runs:-1:1
    if total_runs>=cohort_pars.min_parfor_num
        [idx,out_temp] = fetchNext(F);
    else
        idx = ri;
        [sample_ind,vp_ind{colons{:}}] = ind2sub([nsamps_per_condition,cohort_size],ri);
        if sample_ind == nsamps_per_condition % going in reverse order here, so look for the highest sample ind to reset
            start_ind = 1;
        end
        if cohort.ids(ri)~=""
            out_temp = grabSimData(cohort.ids(ri));
        else
            for vpi = 1:numel(cohort.lattice_parameters)
                M = setField(M,cohort.lattice_parameters(vpi).path,cohort.lattice_parameters(vpi).values(vp_ind{vpi},:));
            end
            [sim_this,sims_to_check,start_ind,sim_id] = findSimilarSims(M,sims_to_check,start_ind);
            if sim_this || any(ids(:,vp_ind{:})==sim_id) % make sure that we didn't already record this sim for these parameter values
                out_temp = simPatient(M);
            else % if the above code did not (somehow) grab this sim, then grab this sim data now
                out_temp = grabSimData(sim_id);
            end
        end
    end
    cohort_pars = updateCohortTimer(cohort_pars,total_runs-ri+1);
    if isfield(out_temp.save_pars,"sim_identifier")
        cohort.ids(idx) = out_temp.save_pars.sim_identifier;
    end
end

if M.save_pars.dt < Inf
    cohort.ids = reshape(cohort.ids,[nsamps_per_condition,cohort_size,1]);
end

cohort.ids = squeeze(permute(cohort.ids,[2:length(cohort_size)+1,1])); % put the sample dimension back along last dimension

if ~isfield(cohort_pars,"cohort_identifier")
    cohort_pars.cohort_identifier = string(datetime("now","Format","yyMMddHHmmssSSS")); % default to this for determining an id if none given
end

while exist(sprintf("data/cohort_%s",cohort_pars.cohort_identifier),"dir") % just in case this directory already exists somehow (not sure how to processes could start at the same time to the millisecond and then one create this folder before the other looks for it)
    cohort_pars.cohort_identifier = string(datetime("now","Format","yyMMddHHmmssSSS")); % default to this for determining an id if none given
end

mkdir(sprintf("data/cohort_%s",cohort_pars.cohort_identifier))

save(sprintf("data/cohort_%s/cohort_%s",cohort_pars.cohort_identifier,cohort_pars.cohort_identifier),"nsamps_per_condition","total_runs","cohort_size")
save(sprintf("data/cohort_%s/cohort_%s",cohort_pars.cohort_identifier,cohort_pars.cohort_identifier),'-struct',"cohort","-append")

end

function lattice_parameters = grabFields(S,lattice_parameters,incoming_struct_path)

fn = string(fieldnames(S));
for i = 1:numel(fn)
    current_struct_path = [incoming_struct_path,fn(i)];
    if isstruct(S.(fn(i)))
        lattice_parameters = grabFields(S.(fn(i)),lattice_parameters,current_struct_path);
    elseif numel(S.(fn(i)))>1 % then vary over these parameters
        lattice_parameters(end+1) = struct("path",current_struct_path,"values",S.(fn(i)));
    end
end

end

function S = setField(S,path,val)

if iscell(path)
    for i = 1:numel(path)
        S = setField(S,path{i},val(i));
    end
elseif length(path)>1
    S.(path(1)) = setField(S.(path(1)),path(2:end),val);
else
    S.(path(1)) = val;
end

end

function [sim_this,sims_to_check,new_start_ind,sim_id] = findSimilarSims(M,sims_to_check,start_ind)
sim_id = ""; 
new_start_ind = numel(sims_to_check)+1; % if no match is found for these pars, then don't search for a match next time
sim_this = true;
for i = start_ind:numel(sims_to_check) % look for the first one that matches the inputs here
    if exist(sprintf("data/%s/output_constants.mat",sims_to_check(i).name),"file") && exist(sprintf("data/%s/output_final.mat",sims_to_check(i).name),"file")
        X = load(sprintf("data/%s/output_constants.mat",sims_to_check(i).name));
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
%                         disp(fn{j})
%                         disp(par_fn{k})
%                         if strcmp(sims_to_check(i).name,"221012070214082")
%                             disp('')
%                         end
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
            sim_id = string(sims_to_check(i).name);
            sims_to_check(i) = [];
            new_start_ind = i;
            return;
        end
    end
end
end

function out = grabSimData(folder_name)

% load(sprintf("data/%s/output_final.mat",folder_name),"tracked")
% out.tracked = tracked;
out.save_pars.sim_identifier = folder_name;

end
