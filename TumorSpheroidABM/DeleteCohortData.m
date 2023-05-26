clearvars;

% This script will fully delete all cohort data, including associated
% simulation data in data/sims. This can be useful after a cohort test run
% whose simulations will not be used subsequently.

delete_data = true; % whether to actually delete these sims, or just identify them (MAKE SURE THEY AREN'T USED BY ANOTHER COHORT!!!!)

cohort_names = "cohort_2305261431";

for j = 1:numel(cohort_names)
    load("data/" + cohort_names(j) + "/output.mat","ids")

    if delete_data
        if input("Are you certain these sims do not belong to another cohort? Cancel deletion with 0 (zero). Proceed with any other numeric input. ")
            for i = 1:numel(ids)
                if exist(sprintf("data/sims/%s",ids(i)),"dir")
                    rmdir(sprintf("data/sims/%s",ids(i)),'s');
                end
            end
        end
    else
        fprintf("Found %d simulations associated with this cohort. They are not deleted though.\n",numel(ids))
    end

    rmdir("data/" + cohort_names(j),'s')
end
