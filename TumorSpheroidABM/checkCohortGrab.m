function checkCohortGrab(cohort_name)

start_wait = tic;
filename = sprintf("data/cohort_%s/output.mat",cohort_name);
while toc(start_wait) < 3 && ~exist(filename,"file")
end
if ~exist(filename,"file")
    error("still waiting for file?")
end
C = load(filename);


C.ids = reshape(C.ids,[],C.nsamps_per_condition);

idx = cell(length(C.cohort_size),1);
for i = 1:size(C.ids,1)

    [idx{:}] = ind2sub(C.cohort_size,i);
    for j = 1:C.nsamps_per_condition
        S = load(sprintf("data/sims/%s/output_constants.mat",C.ids(i,j)));
        for k = 1:numel(C.lattice_parameters)
            if iscell(C.lattice_parameters(k).path)
                for l = 1:numel(C.lattice_parameters(k).path)
                    s = S;
                    for m = 1:length(C.lattice_parameters(k).path{l})
                        s = s.(C.lattice_parameters(k).path{l}(m));
                    end
                    assert(isequal(C.lattice_parameters(k).values(idx{k},l),s))
                end

            else
                s = S;
                for l = 1:length(C.lattice_parameters(k).path)
                    s = s.(C.lattice_parameters(k).path(l));
                end
                assert(isequal(C.lattice_parameters(k).values(idx{k}),s))
            end

        end


    end


end