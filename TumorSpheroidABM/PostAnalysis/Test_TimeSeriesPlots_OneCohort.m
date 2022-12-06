clearvars;
addpath("~/Documents/MATLAB/myfunctions/")

cohort_name = "cohort_660154731186125";
load(sprintf("../data/%s/%s.mat",cohort_name,cohort_name))
nsamps_per_condition = numel(tracked);
f = gobjects(0,1);
%% all tumor sizes
f(end+1,1) = figure;
ax = gca;
hold on
for si = 1:nsamps_per_condition
    plot(tracked(si).t,tracked(si).NT)
end

%% patch tumor size
t = zeros(0,nsamps_per_condition);
NT = zeros(0,nsamps_per_condition);
f(end+1,1) = figure("Name","tumor_size"); hold on;
colors = lines(1);
for si = 1:nsamps_per_condition
    nt = length(tracked(si).t);
    if nt < size(t,1)
        t(1:nt,si) = tracked(si).t;
        t(nt+1:end,si) = NaN;
        NT(1:nt,si) = tracked(si).NT;
        NT(nt+1:end,si) = NaN;
    else
        t(end+1:nt,:) = NaN;
        t(:,si) = tracked(si).t;
        NT(end+1:nt,:) = NaN;
        NT(:,si) = tracked(si).NT;
    end

end

[xx,yy,pc] = my_patchPlot(t,NT,true);
pp = patch(pc{1},max(0,pc{2}),colors,'FaceAlpha',0.15,'EdgeColor','black');
ll = plot(xx,yy,'Color',colors,"LineWidth",1.25);

set(gca,'FontSize',16)
ylabel("Tumor Size (cells)")
xlabel("Time (days)")
f(end).Position = [768   914   632   244];

%% all ctl sizes
f(end+1,1) = figure;
ax = gca;
hold on;

for si = 1:nsamps_per_condition
    plot(tracked(si).t,tracked(si).NI)
end

%% patch ctl size
t = zeros(0,nsamps_per_condition);
NI = zeros(0,nsamps_per_condition);
colors = lines(1);
f(end+1,1) = figure("Name","ctl_size"); hold on;

for si = 1:nsamps_per_condition
    nt = length(tracked(si).t);
    if nt < size(t,1)
        t(1:nt,si) = tracked(si).t;
        t(nt+1:end,si) = NaN;
        NI(1:nt,si) = tracked(si).NI;
        NI(nt+1:end,si) = NaN;
    else
        t(end+1:nt,:) = NaN;
        t(:,si) = tracked(si).t;
        NI(end+1:nt,:) = NaN;
        NI(:,si) = tracked(si).NI;
    end

end

[xx,yy,pc] = my_patchPlot(t,NI,true);
pp = patch(pc{1},max(0,pc{2}),colors,'FaceAlpha',0.15,'EdgeColor','black');
ll = plot(xx,yy,'Color',colors,"LineWidth",1.25);


set(gca,'FontSize',16)
ylabel("CTL Inifiltrate (cells)")
xlabel("Time (days)")
f(end).Position = [768   914   632   244];

%% all ctl proportions

propfn = @(x,y) y./(x+y);
f(end+1,1) = figure;
ax = gca; hold on;
for si = 1:nsamps_per_condition
    plot(tracked(si).t,propfn(tracked(si).NT,tracked(si).NI))
end

%% patch ctl proportions
t = zeros(0,nsamps_per_condition);
Y = zeros(0,nsamps_per_condition);
colors = lines(1);
f(end+1,1) = figure("Name","ctl_proportions"); hold on;

for si = 1:nsamps_per_condition
    nt = length(tracked(si).t);
    if nt < size(t,1)
        t(1:nt,si) = tracked(si).t;
        t(nt+1:end,si) = NaN;
        Y(1:nt,si) = propfn(tracked(si).NT,tracked(si).NI);
        Y(nt+1:end,si) = NaN;
    else
        t(end+1:nt,:) = NaN;
        t(:,si) = tracked(si).t;
        Y(end+1:nt,:) = NaN;
        Y(:,si) = propfn(tracked(si).NT,tracked(si).NI);
    end

end

[xx,yy,pc] = my_patchPlot(t,Y,true);
pp = patch(pc{1},min(1,max(0,pc{2})),colors,'FaceAlpha',0.15,'EdgeColor','none');
ll = plot(xx,yy,'Color',colors,"LineWidth",1.25);

set(gca,'FontSize',16)
ylabel("CTL Proportion")
xlabel("Time (days)")
f(end).Position = [768   914   632   244];

%% print
directories_made = false; % only attempt to make the directories once
file_formats = ["fig","png"];
for i = 1:numel(f)
    if ~isempty(f(i).Name)
        if ~directories_made && ~exist("../figs","dir")
            mkdir("../figs")
        end
        for ffi = 1:length(file_formats)
            if ~directories_made && ~exist(sprintf("../figs/%s/%s",cohort_name,file_formats(ffi)),"dir")
                mkdir(sprintf("../figs/%s/%s",cohort_name,file_formats(ffi)))
            end
            if ~exist(sprintf("../figs/%s/%s/%s.%s",cohort_name,file_formats(ffi),f(i).Name,file_formats(ffi)),"file")
                if file_formats(ffi)=="fig"
                    savefig(f(i),sprintf("../figs/%s/fig/%s",cohort_name,f(i).Name))
                else
                    print(f(i),sprintf("../figs/%s/%s/%s",cohort_name,file_formats(ffi),f(i).Name),sprintf("-d%s",file_formats(ffi)))
                end
            end
        end
        directories_made = true;
    end
end