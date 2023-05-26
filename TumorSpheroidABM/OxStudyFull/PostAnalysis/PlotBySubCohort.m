clearvars;

addpath("../..")
cohort_id = "cohort_2305261503";

load(sprintf("../../data/%s/output.mat",cohort_id),"ids","lattice_parameters")

for ci = 1:3
    figure;
    ax = gobjects(2,2);
    for i = 1:2
        for j = 1:2
            ax(i,j) = subplot(2,2,r2c(2,2,[i,j]),"NextPlot","add");
        end
    end

    cycle = buildCycle();

    for i = 1:numel(ids)
        [i1,i2,i3,~] = ind2sub(size(ids),i);
        if i3~=ci
            continue;
        end
        load(sprintf("../../data/sims/%s/output_final.mat",ids(i)),"tracked")
        load(sprintf("../../data/sims/%s/output_constants.mat",ids(i)))
        plot(ax(i1,i2),tracked.t,tracked.phases(:,cycle.m))
    end
    normalizeYLims(ax)
end

rmpath("../..")
