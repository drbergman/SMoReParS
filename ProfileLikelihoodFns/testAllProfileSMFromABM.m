function f = testAllProfileSMFromABM(profile_file,n_per_fig,sm_par_display_names)

load(profile_file,"out")

npars = size(out,1);
out = reshape(out,npars,[]);
n_abm_vecs = size(out,2);

threshold = chi2inv(0.95,npars);

nfigs = ceil(n_abm_vecs/n_per_fig);
f = gobjects(nfigs,1);
for fi = 1:nfigs
    f(fi) = figure;
    I = (fi-1)*n_per_fig + (1:n_per_fig);
    if fi == nfigs
        I(I>n_abm_vecs) = [];
    end
    ax = gobjects(n_per_fig,npars,nfigs);
    for i = 1:numel(I)
        for j = 1:npars
            if size(out{j,I(i)},1) == 2
                par_ind = 1;
            else
                par_ind = j;
            end
            ax(i,j,fi) = subplot(n_per_fig,npars,r2c(n_per_fig,npars,[i,j]));
            min_val = min(out{j,I(i)}(end,:));
            x = out{j,I(i)}(par_ind,:);
            y = out{j,I(i)}(end,:);
            plot(ax(i,j,fi),x,y);
            yline(ax(i,j,fi),min_val+threshold,"LineStyle","--")
            ax(i,j,fi).YLim(1) = max(ax(i,j,fi).YLim(1),.9*min_val);
            % ax(i,j,fi).YLim(2) = min(max(ax(i,j,fi).YLim(2),min_val+1.1*threshold),min_val+2*threshold);
            ax(i,j,fi).YLim(2) = max(ax(i,j,fi).YLim(2),min_val+1.1*threshold);
        end
    end

    for j = 1:npars
        title(ax(1,j,fi),sm_par_display_names(j))
    end

    for i = 1:numel(I)
        ylabel(ax(i,1,fi),sprintf("#%d",I(i)),"FontWeight","bold")
    end
end
