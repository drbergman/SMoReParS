function suffix = fileSuffix(opts)

suffix = opts.acceptance_method;

if isfield(opts,"acceptance_method") && startsWith(opts.acceptance_method,"specified_parameter_combinations")
    for i = 1:numel(opts.par_combos)
        suffix = suffix + "_";
        for j = 1:numel(opts.par_combos{i})
            suffix = suffix + num2str(opts.par_combos{i}(j));
        end
    end
end

suffix = suffix + ".mat";
