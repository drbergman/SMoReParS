clearvars;

addpath("../..")
cohort_id = "cohort_2305241433";

load(sprintf("../../data/%s/output.mat",cohort_id),"ids","lattice_parameters")

ncohorts = numel(ids)/size(ids,ndims(ids));
x=(1:ncohorts)';
colors = [sin(2*pi/3 * x).^2,cos(2*pi/3 * x).^2,0.5*(sin(2*pi/3 * x).^2+cos(2*pi/3 * x).^2)];

nsamps = numel(ids)/ncohorts;

ax = gobjects(3,1);
for i = 1:3
    ax(i) = subplot(3,1,i,"NextPlot","add");
end

val = buildLatticeVals();

for i = 1:numel(ids)
    [color_ind,~] = ind2sub([ncohorts,nsamps],i);
    load(sprintf("../../data/sims/%s/output_final.mat",ids(i)),"tracked")
    load(sprintf("../../data/sims/%s/output_constants.mat",ids(i)))
    for j = 1:3
        plot(ax(j),tracked.t,tracked.phase_cell_days(:,j) / (tracked.t(2)-tracked.t(1)),'Color',colors(color_ind,:))
        if i==1
            switch j
                case val.phase_g0
                    title(ax(j),"G0")
                case val.phase_g1
                    title(ax(j),"G1")
                case val.phase_m
                    title(ax(j),"M")
            end
        end
    end
end
normalizeYLims(ax)

rmpath("../..")
