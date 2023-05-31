% This script will be a one-off to reduce the size of output_constants.mat
% for each simulation. Instead, all that will be saved there is the
% sim_identifiers the simulations would have based on the cohorts it
% belongs to.

clearvars;

cohorts = dir("data/cohort_*");

for i = 1:numel(cohorts)

    C = load(sprintf("%s/%s/output.mat",cohorts(i).folder,cohorts(i).name),"ids");
    I = strfind(cohorts(i).name,'_');
    cohort_id = cohorts(i).name(I+1:end);

    for j = 1:numel(C.ids)
        file_name = sprintf("data/sims/%s/output_constants.mat",C.ids(j));
        S = load(file_name,"sim_identifier");
        if isfield(S,"sim_identifier")
            sim_identifier = S.sim_identifier;
        else
            sim_identifier = strings(0,1);
        end
        sim_identifier(end+1) = sprintf("%s_%d",cohort_id,j);
        save(file_name,"sim_identifier")

    end


end




