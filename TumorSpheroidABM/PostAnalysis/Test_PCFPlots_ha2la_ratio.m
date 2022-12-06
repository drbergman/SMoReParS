clearvars;
addpath("~/Documents/MATLAB/myfunctions/")

cohort_name = "cohort_69891714061750";
load(sprintf("../data/%s/pcf_low_antigen_to_immune.mat",cohort_name))
out_la = out;
load(sprintf("../data/%s/pcf_high_antigen_to_immune.mat",cohort_name))
out_ha = out;
load(sprintf("../data/%s/%s.mat",cohort_name,cohort_name))
f = gobjects(0,1);

%% all time series plots by assumption
nr = ceil(sqrt(nsamps_per_condition));
nc = ceil(nsamps_per_condition/nr);
ax = gobjects(size(out));
for recruit_ind = 1:2
    for cytotox_ind = 1:2
        f(end+1,1) = figure;
        for si = 1:nsamps_per_condition
            ax(recruit_ind,cytotox_ind,si) = subplot(nr,nc,si);
            pcf = out_ha(recruit_ind,cytotox_ind,si).avg'./out_la(recruit_ind,cytotox_ind,si).avg';
            pcf(out_la(recruit_ind,cytotox_ind,si).avg'==0) = NaN;
            pcfTimeSeriesPlot(options.rr,out(recruit_ind,cytotox_ind,si).t,pcf)
            drawnow
        end
    end
end

%% all time series plots by sample
nr = 2;
nc = 2;
ax = gobjects(size(out));
for si = 1:nsamps_per_condition
    f(end+1) = figure;
    for recruit_ind = 1:2
        for cytotox_ind = 1:2
            ax(recruit_ind,cytotox_ind,si) = subplot(nr,nc,r2c(nr,nc,[recruit_ind,cytotox_ind]));
            pcf = out_ha(recruit_ind,cytotox_ind,si).avg'./out_la(recruit_ind,cytotox_ind,si).avg';
            pcf(out_la(recruit_ind,cytotox_ind,si).avg'==0) = NaN;
            pcfTimeSeriesPlot(options.rr,out(recruit_ind,cytotox_ind,si).t,pcf)
            drawnow
        end
    end
end

%% average time series
nr = 2;
nc = 2;
ax = gobjects(2,2);
f(end+1) = figure("Name","pcf_average_ha2la");
threshold = .0*nsamps_per_condition;
pcf_limits = [Inf,-Inf];
title_strings = ["No effect","Cytotoxic effect";"Recruit effect","Both effects"];
for recruit_ind = 1:2
    for cytotox_ind = 1:2
        ax(recruit_ind,cytotox_ind) = subplot(nr,nc,r2c(nr,nc,[recruit_ind,cytotox_ind]));
        total_avg_ha = zeros(length(options.rr),0);
        total_avg_la = zeros(length(options.rr),0);
        weights = zeros(1,0);
        t = [];
        for si = 1:nsamps_per_condition
            t = unique([t,out(recruit_ind,cytotox_ind,si).t]);
            szs = [size(total_avg_ha,2),size(out(recruit_ind,cytotox_ind,si).avg,2)]; % number of time points for each so far
            if szs(1)<szs(2)
                weights = weights+1;
                total_avg_ha = ((si-1)*total_avg_ha + out_ha(recruit_ind,cytotox_ind,si).avg(:,1:szs(1)))./weights;
                total_avg_la = ((si-1)*total_avg_la + out_la(recruit_ind,cytotox_ind,si).avg(:,1:szs(1)))./weights;
                weights(end+1:szs(2)) = 1;
                total_avg_ha(:,end+1:szs(2)) = out_ha(recruit_ind,cytotox_ind,si).avg(:,szs(1)+1:end);
                total_avg_la(:,end+1:szs(2)) = out_la(recruit_ind,cytotox_ind,si).avg(:,szs(1)+1:end);
            elseif szs(1)==szs(2)
                weights = weights+1;
                total_avg_ha = ((si-1)*total_avg_ha + out_ha(recruit_ind,cytotox_ind,si).avg)./weights;
                total_avg_la = ((si-1)*total_avg_la + out_la(recruit_ind,cytotox_ind,si).avg)./weights;
            else
                weights(1:szs(2)) = weights(1:szs(2))+1;
                total_avg_ha(:,1:szs(2)) = ((si-1)*total_avg_ha(:,1:szs(2)) + out_ha(recruit_ind,cytotox_ind,si).avg)./weights(1:szs(2));
                total_avg_la(:,1:szs(2)) = ((si-1)*total_avg_la(:,1:szs(2)) + out_la(recruit_ind,cytotox_ind,si).avg)./weights(1:szs(2));
            end
        end
        pcf = total_avg_ha(:,weights>threshold)'./total_avg_la(:,weights>threshold)';
        pcf(total_avg_la(:,weights>threshold)'==0) = NaN;
        pcf_limits(1) = min([pcf_limits(1);pcf(:)],[],'omitnan');
        pcf_limits(2) = max([pcf_limits(2);pcf(:)],[],'omitnan');
        pcfTimeSeriesPlot(options.rr,t(weights>threshold),pcf)
        title(title_strings(recruit_ind,cytotox_ind))
    end
end

normalizeYLims(f(end))

pcf_limits(2) = min(pcf_limits(2),4);
one_color = 1*ones(1,3);
ngrid = 1001;
max_colors = 101; % max colors for each of positive and negative directions
positive_corr_n = round(ngrid*(pcf_limits(2)-1)/diff(pcf_limits));
negative_corr_n = ngrid-positive_corr_n;

% unique_pos_colors = floor((max_colors+1)*linspace(0,1,positive_corr_n))/(max_colors+1);
% unique_neg_colors = floor((max_colors+1)*linspace(0,1,negative_corr_n))/(max_colors+1);
unique_pos_color_green = floor((max_colors+1)*linspace(0,one_color(2),positive_corr_n))/(max_colors+1);
unique_pos_color_blue = floor((max_colors+1)*linspace(0,one_color(3),positive_corr_n))/(max_colors+1);
unique_neg_colors_red = floor((max_colors+1)*linspace(0,one_color(1),negative_corr_n))/(max_colors+1);
unique_neg_colors_green = floor((max_colors+1)*linspace(0,one_color(2),negative_corr_n))/(max_colors+1);

% positive_colors = [ones(positive_corr_n,1),flip(unique_pos_colors')*[1,1]];
% negative_colors = [unique_neg_colors'*[1,1],ones(negative_corr_n,1)];
positive_colors = [linspace(1,one_color(1),positive_corr_n)',flip(unique_pos_color_green'),flip(unique_pos_color_blue')];
negative_colors = [unique_neg_colors_red',unique_neg_colors_green',linspace(1,one_color(3),negative_corr_n)'];

cmap = [negative_colors;positive_colors];
for axi = 1:numel(ax)
    colormap(ax(axi),cmap);
    clim(ax(axi),pcf_limits)
    c = colorbar(ax(axi));
    c.Label.String = "g_{ha}(r) / g_{la}(r)";
    c.Label.Rotation = 270;
    c.Label.VerticalAlignment = "baseline";
end
set(ax,'FontSize',16)
f(end).Position = [1231         738         872         420];

% %% plot single sample from each
% nr=2;
% nc=2;
% ax = gobjects(2,2);
% f(end+1) = figure;
% tvals = [2,6,15,25,30,50];
% colors = parula(length(tvals));
% for recruit_ind = 1:2
%     for cytotox_ind = 1:2
%         ax(recruit_ind,cytotox_ind) = subplot(nr,nc,r2c(nr,nc,[recruit_ind,cytotox_ind]));
%         for ti = 1:length(tvals)
%             i = find(out(recruit_ind,cytotox_ind,1).t<=tvals(ti),1,'last');
%             xx = [options.rr,flip(options.rr)]';
%             yy = [out(recruit_ind,cytotox_ind,1).avg(:,i);flip(out(recruit_ind,cytotox_ind,1).avg(:,i))] + [out(recruit_ind,cytotox_ind,1).std(:,i);-flip(out(recruit_ind,cytotox_ind,1).std(:,i))];
%             
%             bad_ind = isnan(xx) | isnan(yy);
%             xx(bad_ind) = [];
%             yy(bad_ind) = [];
%             
%             if isempty(xx)
%                 continue;
%             end
%             pp(recruit_ind,cytotox_ind,ti) = patch(xx,max(0,yy),colors(ti,:),'FaceAlpha',0.2);
%         end
% 
% 
%     end
% end

%% print
directories_made = false; % only attempt to make the directories once
file_formats = ["fig","png"];
for i = 1:numel(f)
    if ishandle(f(i)) && ~isempty(f(i).Name)
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
