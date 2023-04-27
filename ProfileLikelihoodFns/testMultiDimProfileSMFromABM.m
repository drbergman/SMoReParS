function [f,ax,I] = testMultiDimProfileSMFromABM(profile_file,nsamps,input_opts)

opts = defaultTestMultiDimProfileSMFromABMOptions;
if nargin>=3 && ~isempty(input_opts)
    opts = overrideDefaultOptions(opts,input_opts);
end

load(profile_file,"MDProfile","par_names","par_vals","cohort_size","profile_size")

npars = numel(par_names);
n_abm_vecs = prod(cohort_size);

MDProfile = reshape(MDProfile,[n_abm_vecs,profile_size]);
MDProfile = permute(MDProfile,[2:ndims(MDProfile),1]); % vary the ABM par vectors along final dimension

MDProfile = opts.transform_fn(MDProfile);

I = randsample(n_abm_vecs,nsamps,false);

f = figure;
if any(opts.plot_type == ["mesh","contourf"])
    switch opts.plot_type
        case "mesh"
            plotfn = @mesh;
        case "contourf"
            plotfn = @contourf;
    end
    n_sm_par_pairs = nchoosek(npars,2); % number of pairs to do mesh plots
    ax = gobjects(nsamps,n_sm_par_pairs);
    s1 = 1; % index of first sm par to vary (on x-axis)
    s2 = 1; % index of second sm par to vary (on y-axis)
    Idx = cell(npars,1);
    for ci = 1:n_sm_par_pairs % sm parameter to fix for m
        if s2==npars
            s1 = s1+1;
            s2 = s1+1;
        else
            s2 = s2+1;
        end
        x = par_vals{s1};
        y = par_vals{s2};
        for ri = 1:nsamps
            ax(ri,ci) = subplot(nsamps,n_sm_par_pairs,r2c(nsamps,n_sm_par_pairs,[ri,ci]));
            temp = sliceof(MDProfile,ndims(MDProfile),I(ri));
            [~,ind] = max(temp,[],"all");
            [Idx{:}] = ind2sub(profile_size,ind);
            for pi = 1:npars
                if pi~=s1 && pi~=s2
                    temp = sliceof(temp,pi,Idx{pi});
                end
            end
            plotfn(ax(ri,ci),x,y,squeeze(temp)') % transpose because MATLAB uses rows to vary in y direction with mesh
            xlabel(ax(ri,ci),par_names(s1))
            ylabel(ax(ri,ci),par_names(s2))
            if opts.plot_type == "contourf"
                colorbar;
            end
        end
        title(ax(1,ci),sprintf("%s and %s",par_names(s1),par_names(s2)),"FontWeight","bold","FontSize",16)
        if any(par_names(s1)==opts.log_scale_pars)
            set(ax(:,ci),"XScale","log")
        end
        if any(par_names(s2)==opts.log_scale_pars)
            set(ax(:,ci),"YScale","log")
        end
    end
    if opts.plot_type=="mesh"
        set(ax,"ZScale",opts.likelihood_scale)
    elseif opts.plot_type=="contourf"
        set(ax,"ColorScale",opts.likelihood_scale)
    end
elseif opts.plot_type == "lines"
    ax = gobjects(nsamps,npars);
    for ri = 1:nsamps
        S = sliceof(MDProfile,ndims(MDProfile),I(ri));
        for ci = 1:npars
            ax(ri,ci) = subplot(nsamps,npars,r2c(nsamps,npars,[ri,ci]),"NextPlot","add");
            SP = permute(S,[ci,setdiff(1:ndims(S),ci)]);
            SP = reshape(SP,size(SP,1),[]);
            J = randsample(size(SP,2),opts.n_lines,false);
            for j = 1:opts.n_lines
                plot(ax(ri,ci),par_vals{ci},SP(:,J(j)))
            end
        end
    end
    for ci = 1:npars
        title(ax(1,ci),par_names(ci),"FontWeight","bold","FontSize",16)
        if any(par_names(ci)==opts.log_scale_pars)
            set(ax(:,ci),"XScale","log")
        end
    end
    set(ax,"YScale",opts.likelihood_scale)

end

for ri = 1:nsamps
    ylabel(ax(ri,1),sprintf("#%d",I(ri)),"FontWeight","bold")
end


end

function default_options = defaultTestMultiDimProfileSMFromABMOptions

default_options.transform_fn = @(x) x; % whether to use log-likelihood or compute (and use) likelihood
default_options.plot_type = "mesh"; % how to plot; options are currently "mesh" or "lines"
default_options.n_lines = 5; % if using lines, how many lines to plot for each parameter
default_options.likelihood_scale = "linear"; % scale for the log-likelihood axis
default_options.log_scale_pars = []; % list of pars to show in log scale

end