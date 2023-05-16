% A crude way to determine which parameters are most "insensitive"
% This will measure SD when fixing a given parameter at a given value to measure
% the "insensitivity"; combine them to get a total insensitivity score for
% that parameter. The most insensitive parameters will be used to vary
% within subplots so that each subplot should have minimal variability,
% i.e. most of the variation is explained by the parameters defining the
% subplots

clearvars;

save_fig_opts.save_figs = true;
save_fig_opts.file_types = ["fig","png"];
save_fig_opts.fig_names = "ABMOutput";
plot_patches = false;

load("data/summary.mat","D","*par*","t")
Last = sliceof(arrayify(D,"A",1),1,length(t));
Last = squeeze(Last);
insensitivity = cell(4,1);
total_insensitivity = zeros(4,1);
for i = 1:4
    insensitivity{i} = zeros(size(Last,i),1);
    for j = 1:size(Last,i)
        temp = sliceof(Last,i,j);
        insensitivity{i}(j) = std(temp(:));
    end
    total_insensitivity(i) = sum(insensitivity{i}.^2);
end

%% prepare data for subplots
[~,order] = sort(total_insensitivity,1,"descend"); % as similar as possible within a subplot; different as possible across subplot columns (I prefer this to below)
% [~,order] = sort(total_insensitivity,1,"ascend"); % as similar as possible across subplot columns; different as possible within a subplot
D = permute(D,[1;order+1]);
par_names = par_names(order);
display_par_names = display_par_names(order);
par_val_str = par_val_str(order);

%% flip dimension that will vary along rows
D = flip(D,ndims(D)-1);
par_val_str{end-1} = flip(par_val_str{end-1});

%% make subplots
f=figure;

nr = size(D,ndims(D)-1);
nc = size(D,ndims(D));
ax = gobjects(nr,nc);
colors = lines(size(D,3));
lsty = {'-','-.','--',':'};
line_width = 2;
for ri = 1:nr
    for ci = 1:nc
        ax(ri,ci) = subplot(nr,nc,r2c(nr,nc,[ri,ci])); hold on;
        for i = 1:size(D,2)
            for j = 1:size(D,3)
                a = D(1,i,j,ri,ci).A;
                if plot_patches
                    s = D(1,i,j,ri,ci).S; %#ok<UNRCH>
                    patch(ax(ri,ci),[t,flip(t)],[a-s;flip(a+s)],colors(j,:),"EdgeColor","none","FaceAlpha",0.2)
                end
                plot(ax(ri,ci),t,a,"Color",colors(j,:),"LineStyle",lsty{i},"LineWidth",line_width)
            end
        end
    end
end
for ri = 1:nr 
    % set yaxis labels based on the parameter varying across rows
    ylabel(ax(ri,1),sprintf("%s = %s",display_par_names(3),par_val_str{3}(ri)),"FontWeight","bold")
end
for ci = 1:nc
    % set xaxis labels based on the parameter varying across columns
    xlabel(ax(nr,ci),sprintf("%s = %s",display_par_names(4),par_val_str{4}(ci)),"FontWeight","bold")
end
drawnow

% make legend
set(ax,"XLimMode","manual","YLimMode","manual") % do not update axes for these next plots
color_lines = gobjects(size(D,3),1);
for i = 1:size(D,3)
    color_lines(i) = line(ax(1,1),[-2,-1],[0,0],"Color",colors(i,:),"LineWidth",line_width,"DisplayName",sprintf("%s = %s",display_par_names(2),par_val_str{2}(i)));
end
style_lines = gobjects(size(D,2),1);
for i = 1:size(D,2)
    style_lines(i) = line(ax(1,2),[-2,-1],[0,0],"Color",.4*ones(1,3),"LineStyle",lsty{i},"LineWidth",line_width,"DisplayName",sprintf("%s = %s",display_par_names(1),par_val_str{1}(i)));
end

legend(ax(1,1),color_lines,"location","best","FontSize",16)
legend(ax(1,2),style_lines,"location","best","FontSize",16)

% other plot stuff
if t(end)==75
    set(ax,"XTick",0:25:75)
end
set(ax,"FontSize",16)
f.Units = "pixels";
f.Position = [0 0 1474 875];

%% save the figure
saveFigures(f,save_fig_opts)

