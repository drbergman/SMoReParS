clearvars;

cohort_name = "cohort_221217165846308";

C = load(sprintf("data/%s/output.mat",cohort_name));

C.ids = reshape(C.ids,[],C.nsamps_per_condition);

idx = cell(length(C.cohort_size),1);
for i = 1:size(C.ids,1)

    [idx{:}] = ind2sub(C.cohort_size,i);
    for j = 1:C.nsamps_per_condition
        S = load(sprintf("data/sims/%s/output_constants.mat",C.ids(i,j)));
        for k = 1:numel(C.lattice_parameters)
            s = S;
            for l = 1:length(C.lattice_parameters(k).path)
                s = s.(C.lattice_parameters(k).path(l));
            end
            assert(isequal(C.lattice_parameters(k).values(idx{k}),s))


        end


    end


end