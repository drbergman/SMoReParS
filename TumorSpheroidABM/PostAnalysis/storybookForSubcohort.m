function storybookForSubcohort(cohort_name,ids,subcohort_ind,options)

%% set up useful variables
f = gobjects(0,1);
nsamps = numel(ids);
tracked = loadTracked(ids);
t = arrayifyNonuniform(tracked,"t");
t = squeeze(t)';
NT = arrayifyNonuniform(tracked,"NT");
NT = squeeze(NT)';
max_t = ceil(max(t,[],'all','omitnan'));
[r_temp,c_temp] = find(NT~=0);
max_t = min(max_t,ceil(max(t(sub2ind(size(t),r_temp,c_temp)))));

if ~isfield(options,"color")
    options.color = [0 0 0];
end
if ~isfield(options,"type_colors")
    options.type_colors = jet(4);
end
options.type_colors = reshape(options.type_colors,2,2,3);

if isfield(options,"subcohort_title_names")
    name = options.subcohort_title_names(subcohort_ind);
else
    name = string(subcohort_ind);
end

if isfield(options,"subcohort_folder_names")
    folder_name = options.subcohort_folder_names(subcohort_ind);
else
    folder_name = string(subcohort_ind);
end
storybook_folder_path = sprintf("../figs/%s/storybook_%s",cohort_name,folder_name);

if ~isfield(options,"tumor_type_names")
    options.tumor_type_names = ["LA","HA";"LA + Mut","HA + Mut"];
end

if ~isfield(options,"gaussian_kernel_options")
    options.gaussian_kernel_options.h = 4/24; % "width" of the Gaussian
    options.gaussian_kernel_options.n = ceil(max_t*4); % number of time points to sample
end

if ~isfield(options,"reprint")
    options.reprint = false;
end

%% tumor size
fig_name = "tumor_size";
if ~options.reprint && exist(sprintf("%s/fig/%s.fig",storybook_folder_path,fig_name),"file")
    if isfield(options,"redisplay") && options.redisplay
        f(end+1,1) = open(sprintf("%s/fig/%s.fig",storybook_folder_path,fig_name));
    end
else
    f(end+1,1) = figure("Name",fig_name);
    ax = gca;
    ax.NextPlot = "add";
    [xx,yy,pc] = my_patchPlot(t,NT,true);
    pp = patch(ax,pc{1},max(0,pc{2}),options.color,'FaceAlpha',0.15,'EdgeColor','black');
    ll = plot(ax,xx,yy,'Color',options.color,"LineWidth",1.25);
    set(ax,'FontSize',16)
    ylabel(ax,"Tumor Size (cells)")
    xlabel(ax,"Time (days)")
    title(ax,sprintf("%s: Tumor Growth",name))
    xlim(ax,[0 max_t]);
end

%% ctl infiltrate
fig_name = "ctl_infiltrate";
if ~options.reprint && exist(sprintf("%s/fig/%s.fig",storybook_folder_path,fig_name),"file")
    if isfield(options,"redisplay") && options.redisplay
        f(end+1,1) = open(sprintf("%s/fig/%s.fig",storybook_folder_path,fig_name));
    end
else
    f(end+1,1) = figure("Name",fig_name);
    ax = gca;
    ax.NextPlot = "add";
    NI = arrayifyNonuniform(tracked,"NI");
    NI = squeeze(NI)';
    [xx,yy,pc] = my_patchPlot(t,NI,true);
    pp = patch(ax,pc{1},max(0,pc{2}),options.color,'FaceAlpha',0.15,'EdgeColor','black');
    ll = plot(ax,xx,yy,'Color',options.color,"LineWidth",1.25);
    set(ax,'FontSize',16)
    ylabel(ax,"CTL Infiltrate (cells)")
    xlabel(ax,"Time (days)")
    title(ax,sprintf("%s: CTL Infiltrate",name))
    xlim(ax,[0 max_t]);
end

%% ctl proportion
fig_name = "ctl_proportion";
if ~options.reprint && exist(sprintf("%s/fig/%s.fig",storybook_folder_path,fig_name),"file")
    if isfield(options,"redisplay") && options.redisplay
        f(end+1,1) = open(sprintf("%s/fig/%s.fig",storybook_folder_path,fig_name));
    end
else
    f(end+1,1) = figure("Name",fig_name);
    ax = gca;
    ax.NextPlot = "add";
    if ~exist("NT","var")
        NT = arrayifyNonuniform(tracked,"NT");
        NT = squeeze(NT)';
    end
    if ~exist("NIT","var")
        NI = arrayifyNonuniform(tracked,"NI");
        NI = squeeze(NI)';
    end
    P = NI./(NI+NT);
    [xx,yy,pc] = my_patchPlot(t,P,true);
    pp = patch(ax,pc{1},min(1,max(0,pc{2})),options.color,'FaceAlpha',0.15,'EdgeColor','black');
    ll = plot(ax,xx,yy,'Color',options.color,"LineWidth",1.25);
    set(ax,'FontSize',16)
    ylabel(ax,"CTL Proportion")
    xlabel(ax,"Time (days)")
    title(ax,sprintf("%s: CTL Proportion",name))
    xlim(ax,[0 max_t]);
end

%% tumor types
fig_name = "tumor_types";
if ~options.reprint && exist(sprintf("%s/fig/%s.fig",storybook_folder_path,fig_name),"file")
    if isfield(options,"redisplay") && options.redisplay
        f(end+1,1) = open(sprintf("%s/fig/%s.fig",storybook_folder_path,fig_name));
    end
else
    f(end+1,1) = figure("Name",fig_name);
    ax = gca;
    ax.NextPlot = "add";
    pp = gobjects(2,2);
    ll = gobjects(2,2);
    TumTypes = arrayifyNonuniform(tracked,"tumor_types");
    TumTypes = squeeze(TumTypes);
    for mut_ind = 1:2
        for ant_ind = 1:2
            [xx,yy,pc] = my_patchPlot(t,TumTypes(:,:,mut_ind,ant_ind)',true);
            pp(mut_ind,ant_ind) = patch(ax,pc{1},max(0,pc{2}),options.type_colors(mut_ind,ant_ind,:),'FaceAlpha',0.15,'EdgeColor','black',"DisplayName",options.tumor_type_names(mut_ind,ant_ind));
            ll(mut_ind,ant_ind) = plot(ax,xx,yy,'Color',options.type_colors(mut_ind,ant_ind,:),"LineWidth",1.25,"DisplayName",options.tumor_type_names(mut_ind,ant_ind));
        end
    end
    set(ax,'FontSize',16)
    ylabel(ax,"Size (cells)")
    xlabel(ax,"Time (days)")
    title(ax,sprintf("%s: Tumor Size by Type",name))
    xlim(ax,[0 max_t]);
    L = legend(ax,ll(:),"Location","best");
    L.Title.String = "Tumor Type";
end

%% tumor type proportions
fig_name = "tumor_type_proportions";
if ~options.reprint && exist(sprintf("%s/fig/%s.fig",storybook_folder_path,fig_name),"file")
    if isfield(options,"redisplay") && options.redisplay
        f(end+1,1) = open(sprintf("%s/fig/%s.fig",storybook_folder_path,fig_name));
    end
else
    f(end+1,1) = figure("Name",fig_name);
    ax = gca;
    ax.NextPlot = "add";
    pp = gobjects(2,2);
    ll = gobjects(2,2);
    if ~exist("NT","var")
        NT = arrayifyNonuniform(tracked,"NT");
        NT = squeeze(NT)';
    end
    if ~exist("TumTypes","var")
        TumTypes = arrayifyNonuniform(tracked,"tumor_types");
        TumTypes = squeeze(TumTypes);
    end
    for mut_ind = 1:2
        for ant_ind = 1:2
            [xx,yy,pc] = my_patchPlot(t,TumTypes(:,:,mut_ind,ant_ind)'./NT,true);
            pp(mut_ind,ant_ind) = patch(ax,pc{1},min(1,max(0,pc{2})),options.type_colors(mut_ind,ant_ind,:),'FaceAlpha',0.15,'EdgeColor','black',"DisplayName",options.tumor_type_names(mut_ind,ant_ind));
            ll(mut_ind,ant_ind) = plot(ax,xx,yy,'Color',options.type_colors(mut_ind,ant_ind,:),"LineWidth",1.25,"DisplayName",options.tumor_type_names(mut_ind,ant_ind));
        end
    end
    set(ax,'FontSize',16)
    ylabel(ax,"Proportion")
    xlabel(ax,"Time (days)")
    title(ax,sprintf("%s: Tumor Composition",name))
    xlim(ax,[0 max_t]);
    ax.YLim(2) = min(1,ax.YLim(2));
    L = legend(ax,ll(:),"Location","best");
    L.Title.String = "Tumor Type";
end

%% tumor type immune clearance
fig_name = "tumor_type_imm_clearance";
if ~options.reprint && exist(sprintf("%s/fig/%s.fig",storybook_folder_path,fig_name),"file")
    if isfield(options,"redisplay") && options.redisplay
        f(end+1,1) = open(sprintf("%s/fig/%s.fig",storybook_folder_path,fig_name));
    end
else
    f(end+1,1) = figure("Name",fig_name);
    ax = gca;
    ax.NextPlot = "add";
    pp = gobjects(2,2);
    ll = gobjects(2,2);
    ImmClearance = arrayifyNonuniform(tracked,"imm_cleared");
    ImmClearance = squeeze(ImmClearance);
    if ~exist("TumTypes","var")
        TumTypes = arrayifyNonuniform(tracked,"tumor_types");
        TumTypes = squeeze(TumTypes);
    end
    for mut_ind = 1:2
        for ant_ind = 1:2
            for si = nsamps:-1:1
                r(si)=ksr(squeeze(t(1:end-1,si)),ImmClearance(si,2:end,mut_ind,ant_ind)'./(TumTypes(si,1:end-1,mut_ind,ant_ind)'.*diff(t(:,si),1,1)),options.gaussian_kernel_options.h,options.gaussian_kernel_options.n);
            end
            [xx,yy,pc] = my_patchPlot(arrayify(r,"x",1),arrayify(r,"f",1),true);
            pp(mut_ind,ant_ind) = patch(ax,pc{1},max(0,pc{2}),options.type_colors(mut_ind,ant_ind,:),'FaceAlpha',0.15,'EdgeColor','black',"DisplayName",options.tumor_type_names(mut_ind,ant_ind));
            ll(mut_ind,ant_ind) = plot(ax,xx,yy,'Color',options.type_colors(mut_ind,ant_ind,:),"LineWidth",1.25,"DisplayName",options.tumor_type_names(mut_ind,ant_ind));
        end
    end
    set(ax,'FontSize',16)
    ylabel(ax,"Immune Clearance Rate (d^{-1})")
    xlabel(ax,"Time (days)")
    title(ax,sprintf("%s: Immune Clearance by Type",name))
    xlim(ax,[0 max_t]);
    ax.YLim(2) = min(2,ax.YLim(2)); % it seems 2 is a good max value to use otherwise the rate ends up pretty high during the final stages of elimination making it hard to see what happened at early time points
    L = legend(ax,ll(:),"Location","best");
    L.Title.String = "Tumor Type";
end

%% tumor type prolif:contact inhibition
fig_name = "tumor_type_prolif_success_prob";
if ~options.reprint && exist(sprintf("%s/fig/%s.fig",storybook_folder_path,fig_name),"file")
    if isfield(options,"redisplay") && options.redisplay
        f(end+1,1) = open(sprintf("%s/fig/%s.fig",storybook_folder_path,fig_name));
    end
else
    f(end+1,1) = figure("Name",fig_name);
    ax = gca;
    ax.NextPlot = "add";
    pp = gobjects(2,2);
    ll = gobjects(2,2);
    TumProlif = arrayifyNonuniform(tracked,"tum_prolif");
    TumProlif = squeeze(TumProlif);
    TumContactInhibition = arrayifyNonuniform(tracked,"tum_contact_inhibition");
    TumContactInhibition = squeeze(TumContactInhibition);
    Weight = TumProlif + TumContactInhibition;
    Proportion = TumProlif./Weight;
    Avg = sum(Proportion.*Weight,1)./sum(Weight,1);
    Avg = squeeze(Avg);
    STD = zeros(size(Avg));
    for mut_ind = 1:2
        for ant_ind = 1:2
            for ti = 1:size(Proportion,2)
                w = Weight(:,ti,mut_ind,ant_ind);
                w(isnan(w)) = 0;
                STD(ti,mut_ind,ant_ind) = std(Proportion(:,ti,mut_ind,ant_ind),w,'omitnan');
            end
        end
    end
    for mut_ind = 1:2
        for ant_ind = 1:2
            xx = unique(t);
            xx(isnan(xx)) = [];
            xx(end+1:size(t,1)) = NaN;
            yy = Avg(:,mut_ind,ant_ind);
            r_mean = ksr(xx,yy,options.gaussian_kernel_options.h,options.gaussian_kernel_options.n);
            pc{1} = [r_mean.x,flip(r_mean.x)];
            std_temp = STD(:,mut_ind,ant_ind);
            r_std = ksr(xx,std_temp,options.gaussian_kernel_options.h,options.gaussian_kernel_options.n);
            pc{2} = [r_mean.f - r_std.f,flip(r_mean.f+r_std.f)];
            pp(mut_ind,ant_ind) = patch(ax,pc{1},min(1,max(0,pc{2})),options.type_colors(mut_ind,ant_ind,:),'FaceAlpha',0.15,'EdgeColor','black',"DisplayName",options.tumor_type_names(mut_ind,ant_ind));
            ll(mut_ind,ant_ind) = plot(ax,r_mean.x,r_mean.f,'Color',options.type_colors(mut_ind,ant_ind,:),"LineWidth",1.25,"DisplayName",options.tumor_type_names(mut_ind,ant_ind));
        end
    end
    set(ax,'FontSize',16)
    ylabel(ax,"Successful Proliferation Probability")
    xlabel(ax,"Time (days)")
    title(ax,sprintf("%s: Tumor Proliferation Success by Type",name))
    xlim(ax,[0 max_t]);
    L = legend(ax,ll(:),"Location","best");
    L.Title.String = "Tumor Type";
end

%% PCFs:
single_pcf_types = ["","_tumor_to_active","_high_antigen_to_immune","_low_antigen_to_immune","_high_antigen_mut_to_immune","_high_antigen_nonmut_to_immune","_low_antigen_mut_to_immune","_low_antigen_nonmut_to_immune"];
single_title_bits = ["Tumor to Immune","Tumor to Active Immune","HA to Immune","LA to Immune","HA + Mut to Immune","HA + Nonmut to Immune","LA + Mut to Immune","LA + Nonmut to Immune"];
grouped_pcf_types = {["","_tumor_to_active"],["_high_antigen_to_immune","_low_antigen_to_immune"],["_low_antigen_nonmut_to_immune","_low_antigen_mut_to_immune","_high_antigen_nonmut_to_immune","_high_antigen_mut_to_immune"]};
grouped_pcf_filenames = ["tumor_to_immune_types","antigenicity_to_immune","tumor_type_to_immune"];
grouped_title_bits = ["Tumor to Immune Type","Tumor to Immune by Antigenicity","Tumor to Immune by Type"];
grouped_pcf_colors = {lines(2),lines(2),jet(4)};
grouped_legends = {["To Any CTL","To Active CTLs"],["From HA","From LA"],options.tumor_type_names};
for pcf_ind = 1:numel(single_pcf_types)
    temp_name = single_pcf_types(pcf_ind);
    if pcf_ind==1
        temp_name = "_default";
    end
    %% all pcf
    fig_name = sprintf("pcf%s",temp_name);
    if ~options.reprint && exist(sprintf("%s/fig/%s.fig",storybook_folder_path,fig_name),"file")
        if isfield(options,"redisplay") && options.redisplay
            f(end+1,1) = open(sprintf("%s/fig/%s.fig",storybook_folder_path,fig_name));
        end
    elseif ~exist(sprintf("../data/%s/pcf%s.mat",cohort_name,single_pcf_types(pcf_ind)),"file")
        warning("Have not computed some PCF for %s.\n",cohort_name)
    else
        temp = load(sprintf("../data/%s/pcf%s.mat",cohort_name,single_pcf_types(pcf_ind)));
        pcf = temp.out;
        pcf = reshape(pcf,[],nsamps);
        pcf = pcf(subcohort_ind,:);
        pcf_options = temp.options;
        f(end+1,1) = figure("Name",fig_name);
        ax = gca;
        [t_pcf,Avg] = pcfAverage(pcf,pcf_options);
        pcfTimeSeriesPlot(pcf_options.rr,t_pcf,Avg')
        set(ax,'FontSize',16)
        ylabel(ax,"Time (days)")
        xlabel(ax,"r (cell widths)")
        title(ax,sprintf("%s: PCF of %s",name,single_title_bits(pcf_ind)))
    end
end

%% argmax pcf
for group_ind = 1:numel(grouped_pcf_types)
    fig_name = sprintf("pcf_%s_argmax",grouped_pcf_filenames(group_ind));
    if ~options.reprint && exist(sprintf("%s/fig/%s.fig",storybook_folder_path,fig_name),"file")
        if isfield(options,"redisplay") && options.redisplay
            f(end+1,1) = open(sprintf("%s/fig/%s.fig",storybook_folder_path,fig_name));
        end
    else
        clear pcf pcf_options t_pcf argmax rmax; % make sure these variables are not already there
        for pcf_ind = 1:numel(grouped_pcf_types{group_ind})
            skip_this = false;
            if ~exist(sprintf("../data/%s/pcf%s.mat",cohort_name,grouped_pcf_types{group_ind}(pcf_ind)),"file")
                skip_this = true;
                warning("Have not computed some PCF for %s.\n",cohort_name)
                break;
            else
                temp = load(sprintf("../data/%s/pcf%s.mat",cohort_name,grouped_pcf_types{group_ind}(pcf_ind)));
                pcf{pcf_ind} = temp.out;
                pcf{pcf_ind} = reshape(pcf{pcf_ind},[],nsamps);
                pcf{pcf_ind} = pcf{pcf_ind}(subcohort_ind,:);
                pcf_options{pcf_ind} = temp.options;
            end
        end
        if ~skip_this
            f(end+1,1) = figure("Name",fig_name);
            ax = gca;
            ax.NextPlot = "add";
            pp = gobjects(numel(grouped_pcf_types{group_ind}),1);
            ll = gobjects(numel(grouped_pcf_types{group_ind}),1);
            for pcf_ind = 1:numel(grouped_pcf_types{group_ind})
                temp = load(sprintf("../data/%s/pcf%s.mat",cohort_name,grouped_pcf_types{group_ind}(pcf_ind)));
                pcf{pcf_ind} = temp.out;
                pcf{pcf_ind} = reshape(pcf{pcf_ind},[],nsamps);
                pcf{pcf_ind} = pcf{pcf_ind}(subcohort_ind,:);
                pcf_options{pcf_ind} = temp.options;
                t_pcf{pcf_ind} = arrayifyNonuniform(pcf{pcf_ind},"t");
                t_pcf{pcf_ind} = squeeze(t_pcf{pcf_ind});
                for i = 1:numel(pcf{pcf_ind})
                    [~,pcf{pcf_ind}(i).argmax] = max(pcf{pcf_ind}(i).avg,[],1);
                end
                argmax{pcf_ind} = arrayifyNonuniform(pcf{pcf_ind},"argmax");
                argmax{pcf_ind} = squeeze(argmax{pcf_ind});
                rmax{pcf_ind} = NaN(size(argmax{pcf_ind}));
                for i = 1:numel(argmax{pcf_ind})
                    if ~isnan(argmax{pcf_ind}(i))
                        rmax{pcf_ind}(i) = pcf_options{pcf_ind}.rr(argmax{pcf_ind}(i));
                    end
                end
                [xx,yy,pc] = my_patchPlot(t_pcf{pcf_ind}',rmax{pcf_ind}',true);
                pp(pcf_ind) = patch(ax,pc{1},max(0,pc{2}),grouped_pcf_colors{group_ind}(pcf_ind,:),'FaceAlpha',0.15,'EdgeColor','none',"DisplayName",grouped_legends{group_ind}(pcf_ind));
                ll(pcf_ind) = plot(ax,xx,yy,'Color',grouped_pcf_colors{group_ind}(pcf_ind,:),"LineWidth",1.25,"DisplayName",grouped_legends{group_ind}(pcf_ind));
            end
            set(ax,'FontSize',16)
            xlabel(ax,"Time (days)")
            ylabel(ax,"r (cell widths)")
            title(ax,sprintf("%s: Argmax of PCF of %s",name,grouped_title_bits(group_ind)))
            L = legend(ax,ll,"Location","best");
            xlim(ax,[0 max_t])
        end
    end
end

%% rint pcf
for group_ind = 1:numel(grouped_pcf_types)
    fig_name = sprintf("pcf_%s_rint",grouped_pcf_filenames(group_ind));
    if ~options.reprint && exist(sprintf("%s/fig/%s.fig",storybook_folder_path,fig_name),"file")
        if isfield(options,"redisplay") && options.redisplay
            f(end+1,1) = open(sprintf("%s/fig/%s.fig",storybook_folder_path,fig_name));
        end
    else
        clear pcf pcf_options t_pcf rint; % make sure these variables are not already there
        for pcf_ind = 1:numel(grouped_pcf_types{group_ind})
            skip_this = false;
            if ~exist(sprintf("../data/%s/pcf%s.mat",cohort_name,grouped_pcf_types{group_ind}(pcf_ind)),"file")
                skip_this = true;
                warning("Have not computed some PCF for %s.\n",cohort_name)
                break;
            else
                temp = load(sprintf("../data/%s/pcf%s.mat",cohort_name,grouped_pcf_types{group_ind}(pcf_ind)));
                pcf{pcf_ind} = temp.out;
                pcf{pcf_ind} = reshape(pcf{pcf_ind},[],nsamps);
                pcf{pcf_ind} = pcf{pcf_ind}(subcohort_ind,:);
                pcf_options{pcf_ind} = temp.options;
            end
        end
        if ~skip_this
            f(end+1,1) = figure("Name",fig_name);
            ax = gca;
            ax.NextPlot = "add";
            pp = gobjects(numel(grouped_pcf_types{group_ind}),1);
            ll = gobjects(numel(grouped_pcf_types{group_ind}),1);
            for pcf_ind = 1:numel(grouped_pcf_types{group_ind})
                temp = load(sprintf("../data/%s/pcf%s.mat",cohort_name,grouped_pcf_types{group_ind}(pcf_ind)));
                pcf{pcf_ind} = temp.out;
                pcf{pcf_ind} = reshape(pcf{pcf_ind},[],nsamps);
                pcf{pcf_ind} = pcf{pcf_ind}(subcohort_ind,:);
                pcf_options{pcf_ind} = temp.options;
                t_pcf{pcf_ind} = arrayifyNonuniform(pcf{pcf_ind},"t");
                t_pcf{pcf_ind} = squeeze(t_pcf{pcf_ind});
                for i = 1:numel(pcf{pcf_ind})
                    pcf{pcf_ind}(i).rint = sum(pcf{pcf_ind}(i).avg.*(pcf_options{pcf_ind}.rr.^1)',1,'omitnan')./sum(pcf{pcf_ind}(i).avg,1,'omitnan');
                end
                rint{pcf_ind} = arrayifyNonuniform(pcf{pcf_ind},"rint");
                rint{pcf_ind} = squeeze(rint{pcf_ind});
                [xx,yy,pc] = my_patchPlot(t_pcf{pcf_ind}',rint{pcf_ind}',true);
                pp(pcf_ind) = patch(ax,pc{1},max(0,pc{2}),grouped_pcf_colors{group_ind}(pcf_ind,:),'FaceAlpha',0.15,'EdgeColor','none',"DisplayName",grouped_legends{group_ind}(pcf_ind));
                ll(pcf_ind) = plot(ax,xx,yy,'Color',grouped_pcf_colors{group_ind}(pcf_ind,:),"LineWidth",1.25,"DisplayName",grouped_legends{group_ind}(pcf_ind));
            end
            set(ax,'FontSize',16)
            xlabel(ax,"Time (days)")
            ylabel(ax,"r (cell widths)")
            title(ax,sprintf("%s: E[r] of PCF of %s",name,grouped_title_bits(group_ind)))
            L = legend(ax,ll,"Location","best");
            xlim(ax,[0 max_t])
        end
    end
end
% 
%         fig_name = sprintf("pcf%s_rint",temp_name);
%         if ~options.reprint && exist(sprintf("%s/fig/%s.fig",storybook_folder_path,fig_name),"file")
%             if isfield(options,"redisplay") && options.redisplay
%                 f(end+1,1) = open(sprintf("%s/fig/%s.fig",storybook_folder_path,fig_name));
%             end
%         elseif ~exist(sprintf("../data/%s/pcf%s.mat",cohort_name,single_pcf_types(pcf_ind)),"file")
%             warning("Have not computed some PCF for %s.\n",cohort_name)
%         else
%             if ~exist("pcf","var")
%                 temp = load(sprintf("../data/%s/pcf%s.mat",cohort_name,single_pcf_types(pcf_ind)));
%                 pcf = temp.out;
%                 pcf = reshape(pcf,[],nsamps);
%                 pcf = pcf(subcohort_ind,:);
%                 pcf_options = temp.options;
%             end
%             t_pcf = arrayifyNonuniform(pcf,"t");
%             t_pcf = squeeze(t_pcf);
%             for i = 1:numel(pcf)
%                 pcf(i).rint = sum(pcf(i).avg.*(pcf_options.rr.^1)',1,'omitnan')./sum(pcf(i).avg,1,'omitnan');
%             end
%             rint = arrayifyNonuniform(pcf,"rint");
%             rint = squeeze(rint);
%             f(end+1,1) = figure("Name",fig_name);
%             ax = gca;
%             ax.NextPlot = "add";
%             [xx,yy,pc] = my_patchPlot(t_pcf',rint',true);
%             pp = patch(ax,pc{1},max(0,pc{2}),options.color,'FaceAlpha',0.15,'EdgeColor','none');
%             ll = plot(ax,xx,yy,'Color',options.color,"LineWidth",1.25);
%             set(ax,'FontSize',16)
%             xlabel(ax,"Time (days)")
%             ylabel(ax,"r (cell widths)")
%             title(ax,sprintf("%s: E[r] of PCF of %s",name,single_title_bits(pcf_ind)))
%             xlim(ax,[0 max_t])
%         end
% 
%         clear pcf
%     end


% %% PCF: tumor to immune
% fig_name = "pcf_default";
% if ~options.reprint && exist(sprintf("%s/fig/%s.fig",storybook_folder_path,fig_name),"file")
%     if isfield(options,"redisplay") && options.redisplay
%         f(end+1,1) = open(sprintf("%s/fig/%s.fig",storybook_folder_path,fig_name));
%     end
% elseif ~exist(sprintf("../data/%s/pcf.mat",cohort_name),"file") 
%     warning("Have not computed some PCF for %s.\n",cohort_name)
% else
%     temp = load(sprintf("../data/%s/pcf.mat",cohort_name));
%     pcf = temp.out;
%     pcf = reshape(pcf,[],nsamps);
%     pcf = pcf(subcohort_ind,:);
%     pcf_options = temp.options;
%     f(end+1,1) = figure("Name",fig_name);
%     ax = gca;
%     [t_pcf,Avg] = pcfAverage(pcf,pcf_options);
%     pcfTimeSeriesPlot(pcf_options.rr,t_pcf,Avg')
%     set(ax,'FontSize',16)
%     ylabel(ax,"Time (days)")
%     xlabel(ax,"r (cell widths)")
%     title(ax,sprintf("%s: PCF of Tumor to Immune",name))
% end
% 
% %% PCF: tumor to immune arg max
% fig_name = "pcf_default_argmax";
% if ~options.reprint && exist(sprintf("%s/fig/%s.fig",storybook_folder_path,fig_name),"file")
%     if isfield(options,"redisplay") && options.redisplay
%         f(end+1,1) = open(sprintf("%s/fig/%s.fig",storybook_folder_path,fig_name));
%     end
% elseif ~exist(sprintf("../data/%s/pcf.mat",cohort_name),"file") 
%     warning("Have not computed some PCF for %s.\n",cohort_name)
% else
%     if ~exist("pcf","var")
%         temp = load(sprintf("../data/%s/pcf.mat",cohort_name));
%         pcf = temp.out;
%         pcf = reshape(pcf,[],nsamps);
%         pcf = pcf(subcohort_ind,:);
%         pcf_options = temp.options;
%     end
%     t_pcf = arrayifyNonuniform(pcf,"t");
%     t_pcf = squeeze(t_pcf);
%     for i = 1:numel(pcf)
%         [~,pcf(i).argmax] = max(pcf(i).avg,[],1);
%     end
%     argmax = arrayifyNonuniform(pcf,"argmax");
%     argmax = squeeze(argmax);
%     rmax = NaN(size(argmax));
%     for i = 1:numel(argmax)
%         if ~isnan(argmax(i))
%             rmax(i) = pcf_options.rr(argmax(i));
%         end
%     end
%     f(end+1,1) = figure("Name",fig_name);
%     ax = gca;
%     [xx,yy,pc] = my_patchPlot(t_pcf',rmax',true);
%     pp = patch(ax,pc{1},max(0,pc{2}),options.color,'FaceAlpha',0.15,'EdgeColor','none');
%     ll = plot(ax,xx,yy,'Color',options.color,"LineWidth",1.25);
%     set(ax,'FontSize',16)
%     xlabel(ax,"Time (days)")
%     ylabel(ax,"r (cell widths)")
%     title(ax,sprintf("%s: Argmax of PCF of Tumor to Immune",name))
%     xlim(ax,[0 max_t])
% end
% 
% %% PCF: tumor to immune rint
% fig_name = "pcf_default_rint";
% if ~options.reprint && exist(sprintf("%s/fig/%s.fig",storybook_folder_path,fig_name),"file")
%     if isfield(options,"redisplay") && options.redisplay
%         f(end+1,1) = open(sprintf("%s/fig/%s.fig",storybook_folder_path,fig_name));
%     end
% elseif ~exist(sprintf("../data/%s/pcf.mat",cohort_name),"file") 
%     warning("Have not computed some PCF for %s.\n",cohort_name)
% else
%     if ~exist("pcf","var")
%         temp = load(sprintf("../data/%s/pcf.mat",cohort_name));
%         pcf = temp.out;
%         pcf = reshape(pcf,[],nsamps);
%         pcf = pcf(subcohort_ind,:);
%         pcf_options = temp.options;
%     end
%     t_pcf = arrayifyNonuniform(pcf,"t");
%     t_pcf = squeeze(t_pcf);
%     for i = 1:numel(pcf)
%         pcf(i).rint = sum(pcf(i).avg.*(pcf_options.rr.^1)',1,'omitnan')./sum(pcf(i).avg,1,'omitnan');
%     end
%     rint = arrayifyNonuniform(pcf,"rint");
%     rint = squeeze(rint);
%     f(end+1,1) = figure("Name",fig_name);
%     ax = gca;
%     [xx,yy,pc] = my_patchPlot(t_pcf',rint',true);
%     pp = patch(ax,pc{1},max(0,pc{2}),options.color,'FaceAlpha',0.15,'EdgeColor','none');
%     ll = plot(ax,xx,yy,'Color',options.color,"LineWidth",1.25);
%     set(ax,'FontSize',16)
%     xlabel(ax,"Time (days)")
%     ylabel(ax,"r (cell widths)")
%     title(ax,sprintf("%s: E[r] of PCF of Tumor to Immune",name))
%     xlim(ax,[0 max_t])
% end



%% print figures
printFigures(f,cohort_name,storybook_folder_path,options.reprint)
