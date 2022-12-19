function simCohort(M,cohort_pars)

nsamps_per_condition = cohort_pars.nsamps_per_condition;

% lattice sampling of fields of M
fn = string(fieldnames(M));
cohort.lattice_parameters = struct("path",{},"values",{});
for i = 1:numel(fn)
    current_struct_path = fn(i);
    cohort.lattice_parameters = grabFields(M.(fn(i)),cohort.lattice_parameters,current_struct_path);
end

cohort.base_parameters = getBaseParameters(M,cohort.lattice_parameters);


cohort_size = arrayfun(@(i) size(cohort.lattice_parameters(i).values,1),1:numel(cohort.lattice_parameters));
total_runs = prod(cohort_size) * nsamps_per_condition;
cohort_pars.total_runs = total_runs;

colons = repmat({':'},[1,length(cohort_size)]);
vp_ind = cell(1,length(cohort_size));

sims_to_check = dir("data/sims/*");
sims_to_check = sims_to_check([sims_to_check.isdir]);

fn = fieldnames(M);

%% check if previous cohorts ran these sims
previous_cohorts = dir("data/cohort_*");
for i = 1:numel(previous_cohorts)
    no_overlap = true;
    % first check if all non-varied parameters are equal
    PC = load(sprintf("data/%s/output.mat",previous_cohorts(i).name),"ids","lattice_parameters","base_parameters");
    if ~isequal(PC.base_parameters,cohort.base_parameters)
        % then perhaps this is because one varies parameters the other
        % doesn't
        sims_to_check = setdiff(sims_to_check,PC.ids(:));
        continue;
    end
    for vpi = 1:numel(cohort.lattice_parameters)
        for oi = 1:numel(PC.lattice_parameters)
            if isequal(cohort.lattice_parameters(vpi).path,PC.lattice_parameters(oi).path)
                A{vpi} = cohort.lattice_parameters(vpi).values == PC.lattice_parameters(oi).values';
                if all(A{vpi}==false,"all")
                    no_overlap = true;
                end
                break;
            end
        end

        if no_overlap
            break;
        end

    end

    if no_overlap
        continue;
    end

    

end

sims_to_check = sims_to_check(randperm(numel(sims_to_check))); % to not bias samples towards the first sims I ran
cohort.ids = repmat("",[nsamps_per_condition,cohort_size]);

%% record which previous sims have the right parameters
n_found = 0;
for i = 1:numel(sims_to_check)
    if exist(sprintf("data/sims/%s/output_constants.mat",sims_to_check(i).name),"file") && exist(sprintf("data/sims/%s/output_final.mat",sims_to_check(i).name),"file")
        X = load(sprintf("data/sims/%s/output_constants.mat",sims_to_check(i).name));
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

inds_to_run = find(cohort.ids=="");
total_runs = numel(inds_to_run);

%% now fill out the rest of the sim array/grab the sim data identified above
cohort.ids = cohort.ids(:);

if total_runs>=cohort_pars.min_parfor_num
    F(1:total_runs) = parallel.FevalFuture;
    ppool = gcp;
    cohort_pars.num_workers = ppool.NumWorkers;
    for ri = 1:total_runs % run index
        fprintf("Setting up simulation %d of %d...\n",ri,total_runs)
        [~,vp_ind{colons{:}}] = ind2sub([nsamps_per_condition,cohort_size],inds_to_run(ri));
        for vpi = 1:numel(cohort.lattice_parameters)
            M = setField(M,cohort.lattice_parameters(vpi).path,cohort.lattice_parameters(vpi).values(vp_ind{vpi},:));
        end
        M.save_pars.idx_in_cohort = ri;
        F(ri) = parfeval(ppool,@simPatient,1,M);
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
        [temp,out_temp] = fetchNext(F);
        idx = inds_to_run(temp);
    else
        idx = inds_to_run(ri);
        [~,vp_ind{colons{:}}] = ind2sub([nsamps_per_condition,cohort_size],idx);
        for vpi = 1:numel(cohort.lattice_parameters)
            M = setField(M,cohort.lattice_parameters(vpi).path,cohort.lattice_parameters(vpi).values(vp_ind{vpi},:));
        end
        M.save_pars.idx_in_cohort = ri;
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

cohort.ids = squeeze(permute(cohort.ids,[2:length(cohort_size)+1,1])); % put the sample dimension back along last dimension

if ~isfield(cohort_pars,"cohort_identifier")
    cohort_pars.cohort_identifier = string(datetime("now","Format","yyMMddHHmmssSSS")); % default to this for determining an id if none given
end

while exist(sprintf("data/cohort_%s",cohort_pars.cohort_identifier),"dir") % just in case this directory already exists somehow (not sure how to processes could start at the same time to the millisecond and then one create this folder before the other looks for it)
    cohort_pars.cohort_identifier = string(datetime("now","Format","yyMMddHHmmssSSS")); % default to this for determining an id if none given
end

mkdir(sprintf("data/cohort_%s",cohort_pars.cohort_identifier))

save(sprintf("data/cohort_%s/output",cohort_pars.cohort_identifier),"nsamps_per_condition","total_runs","cohort_size")
save(sprintf("data/cohort_%s/output",cohort_pars.cohort_identifier),'-struct',"cohort","-append")

fprintf("Finished cohort. Folder is: cohort_%s\n",cohort_pars.cohort_identifier)
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

function base_parameters = getBaseParameters(M,lattice_parameters)

base_parameters = M;
for i = 1:numel(lattice_parameters)
    base_parameters = prune(base_parameters,lattice_parameters(i).path);
end
end

function s = prune(s,path)
if length(path)>1
    s.(path(1)) = prune(s.(path(1)),path(2:end));
else
    s = rmfield(s,path(1));
end
end
