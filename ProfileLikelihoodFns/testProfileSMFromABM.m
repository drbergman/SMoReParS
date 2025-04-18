function [f,ax,I] = testProfileSMFromABM(profile_file,nsamps,sm_par_display_names,input_opts)

opts = defaultTestProfileSMFromABMOptions;
if nargin>3 && ~isempty(input_opts)
    opts = overrideDefaultOptions(opts,input_opts);
end

load(profile_file,"profiles")

npars = size(profiles,1);
profiles = reshape(profiles,npars,[]);
n_abm_vecs = size(profiles,2);

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
threshold = chi2inv(0.95,npars);
f = figure;
if nsamps>1 || isempty(opts.parameter_layout)
    ax = gobjects(nsamps,npars);
    for i = 1:nsamps
        for j = 1:npars
            if size(profiles{j,I(i)},1) == 2
                par_ind = 1;
            else
                par_ind = j;
            end
            ax(i,j) = subplot(nsamps,npars,r2c(nsamps,npars,[i,j]));
            min_val = min(profiles{j,I(i)}(end,:));
            x = profiles{j,I(i)}(par_ind,:);
            y = profiles{j,I(i)}(end,:);
            plot(ax(i,j),x,y,"LineWidth",opts.LineWidth,"Color",opts.LineColor);
            yline(ax(i,j),min_val+threshold,"LineStyle","--","LineWidth",opts.LineWidth)
            ax(i,j).YLim(1) = max(ax(i,j).YLim(1),.9*min_val);
            ax(i,j).YLim(2) = min(max(ax(i,j).YLim(2),min_val+1.1*threshold),min_val+2*threshold);
        end
    end

    for j = 1:npars
        switch opts.place_par_names
            case "title"
                title(ax(1,j),sm_par_display_names(j),"FontWeight","bold","FontSize",16)
            case "xlabel"
                xlabel(ax(end,j),sm_par_display_names(j),"FontWeight","bold","FontSize",16)
        end
    end

    if opts.show_y_label
        for i = 1:nsamps
            ylabel(ax(i,1),sprintf("#%04d",I(i)),"FontWeight","bold")
        end
    end

else % one sample (so likely a profile from data) and a given parameter layout
    nr = opts.parameter_layout(1);
    nc = opts.parameter_layout(2);
    ax = gobjects(nr,nc);
    for j = 1:npars
        if size(profiles{j,I},1) == 2
            par_ind = 1;
        else
            par_ind = j;
        end
        ax(j) = subplot(nr,nc,j);
        min_val = min(profiles{j,I}(end,:));
        x = profiles{j,I}(par_ind,:);
        y = profiles{j,I}(end,:);
        plot(ax(j),x,y,"LineWidth",opts.LineWidth,"Color",opts.LineColor);
        yline(ax(j),min_val+threshold,"LineStyle","--","LineWidth",opts.LineWidth)
        ax(j).YLim(1) = max(ax(j).YLim(1),.9*min_val);
        ax(j).YLim(2) = min(max(ax(j).YLim(2),min_val+1.1*threshold),min_val+2*threshold);
        switch opts.place_par_names
            case "title"
                title(ax(j),sm_par_display_names(j),"FontWeight","bold","FontSize",16)
            case "xlabel"
                xlabel(ax(j),sm_par_display_names(j),"FontWeight","bold","FontSize",16)
        end
    end
    ax = ax'; % to match the layout of the subplots

end
end

function default_options = defaultTestProfileSMFromABMOptions

default_options.abm_vec_inds = []; % sample from these ABM vecs (empty ==> sample from all)
default_options.LineWidth = 2;
default_options.place_par_names = "title";
default_options.LineColor = lines(1);

default_options.parameter_layout = []; % rows x columns of parameter layout

default_options.show_y_label = true;

end
