clearvars;

addpath("../..")
addpath("~/Documents/MATLAB/myfunctions/")
cohort_id = "cohort_230124175743017";

load(sprintf("../../data/%s/output.mat",cohort_id),"ids","lattice_parameters","nsamps_per_condition")

ncohorts = numel(ids)/nsamps_per_condition;
npars = numel(lattice_parameters);
x=(1:ncohorts)';
val_color = lines(max(size(ids,1:npars)));

figure;
ax = gobjects(npars,5);
for i = 1:5
    for j = 1:npars
        ax(j,i) = subplot(npars,5,r2c(npars,5,[j,i]),"NextPlot","add");
    end
end

cycle = buildCycle();
S_max = 5;

for j = 1:npars
    l_temp = gobjects(size(ids,j),S_max,4);
    for k = 1:size(ids,j)
        id_temp = sliceof(ids,j,k); % grab the ids for sims with the kth value of the jth parameter
        id_temp = id_temp(:);
        S = min(S_max,numel(id_temp));
        id_temp = id_temp(randperm(numel(id_temp),S)); % sample S of these
        for l = 1:S

            [color_ind,~] = ind2sub([ncohorts,nsamps_per_condition],i);
            load(sprintf("../../data/sims/%s/output_final.mat",id_temp(l)),"tracked")
            load(sprintf("../../data/sims/%s/output_constants.mat",id_temp(l)))
            n_phases = size(tracked.phases,2);
            for m = 1:n_phases
                l_temp(k,l,m) = plot(ax(j,m),tracked.t,tracked.phases(:,m),'Color',[val_color(k,:),0.2]);
                l_temp(k,l,m).DisplayName = sprintf("%3.2f",lattice_parameters(j).values(k));
                if l==1
                    if j==1
                        switch m
                            case cycle.g1
                                title(ax(j,m),"G1")
                            case cycle.s
                                title(ax(j,m),"S")
                            case cycle.g2
                                title(ax(j,m),"G2")
                            case cycle.m
                                title(ax(j,m),"M")
                        end
                    end
                    if m==1 && k==size(ids,j)
                        if iscell(lattice_parameters(j).path)
                            ylabel(ax(j,m),regexprep(lattice_parameters(j).path{1}(end),"_"," "))
                        else
                            ylabel(ax(j,m),regexprep(lattice_parameters(j).path(end),"_"," "))
                        end
                        legend(ax(j,m),l_temp(:,1,1),"Location","northwest","AutoUpdate","off")
                    end
                end
            end
            plot(ax(j,n_phases+1),tracked.t,sum(tracked.phases,2),'Color',[val_color(k,:),0.2]);
            if l==1 && j==1
                title(ax(j,n_phases+1),"Total")
            end
        end
    end
end

for m = 1:n_phases
normalizeYLims(ax(:,m))
end

rmpath("../..")
