function [f,I] = testProfileSMFromABM(profile_file,nsamps,sm_par_display_names,input_opts)

opts = defaultTestProfileSMFromABMOptions;
if nargin>3 && ~isempty(input_opts)
    opts = overrideDefaultOptions(opts,input_opts);
end

load(profile_file,"out")

npars = size(out,1);
out = reshape(out,npars,[]);
n_abm_vecs = size(out,2);

if isempty(opts.abm_vec_inds)
    I = randsample(n_abm_vecs,nsamps,false);
else
    opts.abm_vec_inds(opts.abm_vec_inds>n_abm_vecs) = []; % do not keep any indices that exceed the number of abm parameter vectors
    nsamps = min(nsamps,numel(opts.abm_vec_inds)); % limit the number of samples to those selected in options
    I = randsample(opts.abm_vec_inds,nsamps,false);
end
I = unique(I);
nsamps = numel(I);
threshold = chi2inv(0.95,npars);
f = figure;
ax = gobjects(nsamps,npars);
for i = 1:nsamps
    for j = 1:npars
        if size(out{j,I(i)},1) == 2
            par_ind = 1;
        else
            par_ind = j;
        end
        ax(i,j) = subplot(nsamps,npars,r2c(nsamps,npars,[i,j]));
        min_val = min(out{j,I(i)}(end,:));
        x = out{j,I(i)}(par_ind,:);
        y = out{j,I(i)}(end,:);
        plot(ax(i,j),x,y,"LineWidth",2);
        yline(ax(i,j),min_val+threshold,"LineStyle","--")
        ax(i,j).YLim(1) = max(ax(i,j).YLim(1),.9*min_val);
        ax(i,j).YLim(2) = min(max(ax(i,j).YLim(2),min_val+1.1*threshold),min_val+2*threshold);
    end
end

for j = 1:npars
    title(ax(1,j),sm_par_display_names(j),"FontWeight","bold","FontSize",16)
end

for i = 1:nsamps
    ylabel(ax(i,1),sprintf("#%d",I(i)),"FontWeight","bold")
end

end

function default_options = defaultTestProfileSMFromABMOptions

default_options.abm_vec_inds = []; % sample from these ABM vecs (empty ==> sample from all)

end
