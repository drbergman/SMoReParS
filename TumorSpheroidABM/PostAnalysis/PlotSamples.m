clearvars;

addpath("..")
addpath("~/Documents/MATLAB/myfunctions/")
cohort_id = "cohort_221215122756612";

load(sprintf("../data/%s/output.mat",cohort_id),"ids","lattice_parameters")

ncohorts = numel(ids)/size(ids,ndims(ids));
npars = ndims(ids)-1;
x=(1:ncohorts)';
val_color = lines(max(size(ids,1:ndims(ids)-1)));

nsamps = numel(ids)/ncohorts;

figure;
ax = gobjects(npars,3);
for i = 1:3
    for j = 1:npars
        ax(j,i) = subplot(npars,3,r2c(npars,3,[j,i]),"NextPlot","add");
    end
end

val = buildLatticeVals();
S = 100;

for j = 1:npars
    l_temp = gobjects(size(ids,j),S,3);
    for k = 1:size(ids,j)
        id_temp = sliceof(ids,j,k); % grab the ids for sims with the kth value of the jth parameter
        id_temp = id_temp(:);
        id_temp = id_temp(randperm(numel(id_temp),S)); % sample S of these
        for l = 1:S

            [color_ind,~] = ind2sub([ncohorts,nsamps],i);
            load(sprintf("../data/%s/output_final.mat",id_temp(l)),"tracked")
            load(sprintf("../data/%s/output_constants.mat",id_temp(l)))
            for m = 1:3
                l_temp(k,l,m) = plot(ax(j,m),tracked.t,tracked.phase_cell_days(:,m) / (tracked.t(2)-tracked.t(1)),'Color',[val_color(k,:),0.2]);
                l_temp(k,l,m).DisplayName = sprintf("%3.2f",lattice_parameters(j).values(k));
                if l==1
                    if j==1
                        switch m
                            case val.phase_g0
                                title(ax(j,m),"G0")
                            case val.phase_g1
                                title(ax(j,m),"G1")
                            case val.phase_m
                                title(ax(j,m),"M")
                        end
                    end
                    if m==1 && k==size(ids,j)
                        ylabel(ax(j,m),regexprep(lattice_parameters(j).path(end),"_"," "))
                        legend(ax(j,m),l_temp(:,1,1),"Location","northwest","AutoUpdate","off")
                    end
                end
            end
        end
    end
end

for m = 1:3
normalizeYLims(ax(:,m))
end
