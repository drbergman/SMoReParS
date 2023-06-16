function [f,I] = testSMFitToABM(files,nsamps,fn,fn_opts,input_opts)


opts = defaultTestSMFitToABMOptions;
if nargin > 4 && ~isempty(input_opts)
    opts = overrideDefaultOptions(opts,input_opts);
end

if contains(path,"myfunctions")
    path_changed = false;
else
    path_changed = true;
    addpath("~/Documents/MATLAB/myfunctions/")
end

if isfield(files,"optimal_parameters")
    load(files.optimal_parameters,"P")
else
    warning("Rename this field in files from par_file --> optimal_parameters")
    load(files.par_file,"P")
end
if isfield(files,"data")
    load(files.data,"t","D","C","cohort_size","nsamps_per_parameter_vector","n_time_series","n_conditions"); %#ok<NASGU> % some of these variables are not used now, but they might be once I get to filling out the conditional statements below
else
    warning("Rename this field in files from data_file --> data")
    load(files.data_file,"t","D","C","cohort_size","nsamps_per_parameter_vector","n_time_series","n_conditions"); %#ok<NASGU> % some of these variables are not used now, but they might be once I get to filling out the conditional statements below
end

n_sm_pars = size(P,1);
P = reshape(P,n_sm_pars,[]);
n_abm_vecs = size(P,2);

if isempty(opts.abm_vec_inds)
    nsamps = min(nsamps,n_abm_vecs);
    I = randsample(n_abm_vecs,nsamps,false);
else
    opts.abm_vec_inds(opts.abm_vec_inds>n_abm_vecs) = []; % do not keep any indices that exceed the number of abm parameter vectors
    if nsamps < numel(opts.abm_vec_inds)
        I = randsample(opts.abm_vec_inds,nsamps,false);
    elseif nsamps == numel(opts.abm_vec_inds)
        I = opts.abm_vec_inds;
    else
        remaining_population = setdiff(1:n_abm_vecs,opts.abm_vec_inds);
        rand_I = randsample(remaining_population,nsamps-numel(opts.abm_vec_inds),false);
        I = [opts.abm_vec_inds(:);rand_I(:)];
    end
end

I = unique(I);
nsamps = numel(I);
D = reshape(D,n_conditions,[]);

f = gobjects(2,1);
%% these will be filled out as new cases are covered with this
if n_time_series == 1 && n_conditions == 1 % just plot in a rough square
    f(1)=figure;
    nr = ceil(sqrt(nsamps));
    nc = ceil(nsamps/nr);
    ax = gobjects(nsamps,1);
    for i = 1:nsamps
        ax(i) = subplot(nr,nc,i,"NextPlot","add");
        patch(ax(i),[t(:);flip(t(:))],[D(I(i)).A-D(I(i)).S;flipud(D(I(i)).A+D(I(i)).S)],opts.data_color,"FaceAlpha",0.2,"EdgeColor","none")
        plot(ax(i),t,D(I(i)).A,opts.data_color,"LineStyle","--")
        out = fn(P(:,I(i)),t,C{1},fn_opts);
        plot(ax(i),t,out,"-","LineWidth",2,"Color",opts.fit_color)
        title(ax(i),sprintf("#%04d",I(i)),"FontWeight","bold")
    end
elseif n_time_series > 1 && n_conditions == 1
    if isempty(opts.column_names)
        opts.column_names = cell(n_time_series,1);
        for i = 1:n_time_series
            opts.column_names{i} = sprintf("Time Series #%d",i);
        end
    end
    f(1)=figure;
    nr = nsamps;
    nc = n_time_series;
    ax = gobjects(nsamps,n_time_series);
    for ri = 1:nsamps
        out = fn(P(:,I(ri)),t,C{1},fn_opts,D(I(ri)).A);
        for ci = 1:n_time_series
            ax(ri,ci) = subplot(nr,nc,r2c(nr,nc,[ri,ci]),"NextPlot","add");
            patch(ax(ri,ci),[t(:);flip(t(:))],[D(I(ri)).A(:,ci)-D(I(ri)).S(:,ci);flipud(D(I(ri)).A(:,ci)+D(I(ri)).S(:,ci))],opts.data_color,"FaceAlpha",0.2,"EdgeColor","none")
            plot(ax(ri,ci),t,D(I(ri)).A(:,ci),opts.data_color,"LineStyle","--")
            plot(ax(ri,ci),t,out(:,ci),"-","LineWidth",2,"Color",opts.fit_color)
        end
        ylabel(ax(ri,1),sprintf("#%04d",I(ri)),"FontWeight","bold")
    end
    for ci = 1:n_time_series
        title(ax(1,ci),opts.column_names{ci})
    end
elseif n_time_series==2 && n_conditions > 1
    f(1)=figure;
    nr = nsamps;
    nc = n_conditions*n_time_series;
    ax = gobjects(nsamps,n_conditions,n_time_series);

    if isempty(opts.column_names)
        opts.column_names = cell(n_conditions,n_time_series);
        for ci = 1:n_conditions
            for tsi = 1:n_time_series
                opts.column_names{ci,tsi} = sprintf("Condition #%d, Series #%d",ci,tsi);
            end
        end
    end

    for ri = 1:nsamps
        for ci = 1:n_conditions
            out = fn(P(:,I(ri)),t,C{ci},fn_opts);
            for tsi = 1:n_time_series
                ax(ri,ci,tsi) = subplot(nr,nc,r2c(nr,nc,[ri,(ci-1)*n_time_series+tsi]),"NextPlot","add");
                patch(ax(ri,ci,tsi),[t(:);flip(t(:))],[D(ci,I(ri)).A(:,tsi)-D(ci,I(ri)).S(:,tsi);flipud(D(ci,I(ri)).A(:,tsi)+D(ci,I(ri)).S(:,tsi))],opts.data_color,"FaceAlpha",0.2,"EdgeColor","none")
                plot(ax(ri,ci,tsi),t,D(ci,I(ri)).A(:,tsi),opts.data_color,"LineStyle","--")
                plot(ax(ri,ci,tsi),t,out(:,tsi),"-","LineWidth",2,"Color",opts.fit_color)
            end
        end
        ylabel(ax(ri,1,1),sprintf("#%04d",I(ri)),"FontWeight","bold")
    end
    for ci = 1:n_conditions
        for tsi = 1:n_time_series
            title(ax(1,ci,tsi),opts.column_names{ci,tsi})
        end
    end
end
xlim(ax,[t(1) t(end)])
%% histograms of parameter values
f(2)=figure;
nr = ceil(sqrt(n_sm_pars));
nc = ceil(n_sm_pars/nr);
for i = 1:n_sm_pars
    subplot(nr,nc,i);
    hold on
    xL = quantile(P(i,:),[.025,.975]);
    xmin = xL(1);
    xmax = xL(2);
    histogram(P(i,:))
    yL = ylim;
    patch([xmin,xmin,xmax,xmax],[yL,flip(yL)],[.1 .1 .1],"FaceAlpha",.2,"EdgeColor","none")
    switch opts.place_par_names
        case "title"
            title(opts.par_names(i))
        case "xlabel"
            xlabel(opts.par_names(i))
    end
    xline(median(P(i,:)),"LineWidth",2)
end

%% histograms of RSS values
if isfield(files,"sm_fit_file")
    warning("off",'MATLAB:load:variableNotFound')
    load(files.sm_fit_file,"raw_error_options","weights")
    warning("on",'MATLAB:load:variableNotFound')
    if ~exist("raw_error_options","var")
        raw_error_options = [];
    end
    if ~exist("weights","var")
        weights = ones(n_conditions,1);
    end
    f(3) = figure;
    RSS = zeros(1,size(P,2));
    for i = 1:size(P,2)
        RSS(i) = arrayfun(@(j) rawError(P(:,i),t,D(j,i),fn,C{j},fn_opts,raw_error_options),1:n_conditions)*weights;
    end
    if opts.rss_on_log_scale
        RSS = log(RSS);
    end
    histogram(RSS,"Orientation",opts.rss_orientation,"Normalization",opts.rss_normalization);
    if opts.show_rss_smoothing
        mu = mean(RSS);
        sigma = std(RSS);
        [~,edges] = histcounts(RSS);
        hold on
        xx = linspace(min(RSS),max(RSS),1001);
        yy = normpdf(xx,mu,sigma);
        if strcmpi(opts.rss_orientation,"vertical")
            plot(xx,yy*length(RSS)*(edges(2)-edges(1)),"LineWidth",2)
            rss_lim = xlim;
        else
            plot(yy*length(RSS)*(edges(2)-edges(1)),xx,"LineWidth",2)
            rss_lim = ylim;
        end
        if opts.rss_on_log_scale
            rsslog10L = rss_lim*log10(exp(1));
            rss_ticks = ceil(rsslog10L(1)):floor(rsslog10L(2));
            rss_tick_labels = strings(length(rss_ticks),1);
            for i = 1:length(rss_ticks)
                rss_tick_labels(i) = sprintf("10^{%d}",rss_ticks(i));
            end
        end
    end
    if opts.rss_on_log_scale
        if strcmpi(opts.rss_orientation,"vertical")
            xticks(rss_ticks*log(10));
            xticklabels(rss_tick_labels)
            xlabel("RSS")
            ylabel("Frequency")
        else
            yticks(rss_ticks*log(10));
            yticklabels(rss_tick_labels)
            ylabel("RSS")
            xlabel("Frequency")
        end
    end
end


if path_changed
    rmpath("~/Documents/MATLAB/myfunctions/")
end

end

function default_options = defaultTestSMFitToABMOptions

default_options.par_names = []; % 
default_options.column_names = [];
default_options.place_par_names = "title";

default_options.abm_vec_inds = []; % abm vec inds to force into the sample

default_options.rss_on_log_scale = false;
default_options.rss_orientation = "Vertical";
default_options.rss_normalization = "count";
default_options.show_rss_smoothing = false;

%% colors
default_options.data_color = "black";
default_options.fit_color = "red";



end
