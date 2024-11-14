clearvars;

addpath("~/Documents/MATLAB/myfunctions/") % replace with path (rel or abs) to myfunctions
cohort_name = "cohort_221217222613078";
C = load(sprintf("../data/%s/output.mat",cohort_name));

%%
for i = numel(C.ids):-1:1
    S = load(sprintf("../data/sims/%s/output_final.mat",C.ids(i)));
    count(:,i) = S.tracked.NT;
    phase_count(:,:,i) = S.tracked.phases;
end

t = S.tracked.t;
nt = length(t);

%%
count = reshape(count,[nt,size(C.ids)]);
phase_count = reshape(phase_count,[nt,4,size(C.ids)]);

%%
average_count = mean(count,ndims(count));
phase_average = mean(phase_count,ndims(phase_count));

%%
count_std = std(count,[],ndims(count));
phase_std = std(phase_count,[],ndims(phase_count));

%%
% save(sprintf("../data/%s/summary.mat",cohort_name),"average_count","phase_average","count_std","phase_std","-v7.3")

%%
average_count = reshape(average_count,nt,3,3,[]);
count_std = reshape(count_std,nt,3,3,[]);

figure;
for i = 1:3
    for j = 1:3
        subplot(3,3,r2c(3,3,[i,j]),"NextPlot","add")
        for k = 1:size(average_count,4)
            patch([t;flip(t)],[average_count(:,i,j,k)-count_std(:,i,j,k);flip(average_count(:,i,j,k)+count_std(:,i,j,k))],rand(1,3),"FaceAlpha",0.15,"EdgeColor","none");
        end
    end
end

%%
phase_average = reshape(phase_average,nt,4,3,3,[]);
phase_std = reshape(phase_std,nt,4,3,3,[]);

for pj = 1:4
    figure;
    for i = 1:3
        for j = 1:3
            subplot(3,3,r2c(3,3,[i,j]),"NextPlot","add")
            for k = 1:size(phase_average,ndims(phase_average))
                patch([t;flip(t)],[phase_average(:,pj,i,j,k)-phase_std(:,pj,i,j,k);flip(phase_average(:,pj,i,j,k)+phase_std(:,pj,i,j,k))],rand(1,3),"FaceAlpha",0.15,"EdgeColor","none");
            end
        end
    end
end
