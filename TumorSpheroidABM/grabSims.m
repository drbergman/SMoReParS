function [cohort,sims_to_check,finished] = grabSims(cohort,PC,sims_to_check)

finished = false;
all_fn = fieldnames(cohort.all_parameters);
if ~isequal(sort(string(all_fn)),sort(string(fieldnames(PC.all_parameters))))
    % these did not have the same parameters coming in, so move on
    return;
end

[cohort,sims_to_check,finished,~] = continueOn(cohort,PC,sims_to_check,rem_names);
end

function [cohort,sims_to_check,skip_this_cohort] = continueOn(cohort,PC,sims_to_check,rem_names)

skip_this_cohort = false;
if ~isempty(rem_names)
    current_val = aps.(rem_names{1});
    previous_val = PC.all_parameters.(rem_names{1});

    if size(current_val,1)==1 % then this parameter is not currently being varied, check that at least one of previous used it
        if size(previous_val,1)==1 % then this was also not varied previously
            if ~isequal(current_val,previous_val)
                skip_this_cohort = true;
                sims_to_check = 
                return;
            end
        else % then need this was varied previously, need to see if the current val matches any of the previous, and then record where those are
            [used_prev,prev_ind] = ismember(current_val,previous_val,"rows");
            if ~used_prev % make sure the current value was one of the previously used
                skip_this_cohort = true;
                return;
            end
            for k = 1:numel(PC.lattice_parameters)
                if strcmp(rem_names{j},[PC.lattice_parameters(k).path{1},'_DOT_',PC.lattice_parameters(k).path{2}])
                    prev_dim = k;
                    break;
                end
            end
            sims_to_check = setdiff(sims_to_check,removeslice(PC.ids,prev_dim,prev_ind));
            target_dims{prev_dim} = k;
            %                 PC.ids = sliceof(PC.ids,prev_dim,prev_ind);
        end
    end

else % then we should just be at the samples

end
end