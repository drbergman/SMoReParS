function [f,I] = testProfileSMFromABM(profile_file,nsamps,sm_par_display_names,is_cleaned)

load(profile_file,"out")

npars = size(out,1);
out = reshape(out,npars,[]);
n_abm_vecs = size(out,2);

I = randsample(n_abm_vecs,nsamps,false);
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
        if ~is_cleaned
            [par_min,par_max] = getProfileBounds([x;y],threshold);
            inds = x>=par_min & x<=par_max;
            x = x(inds);
            y = y(inds);
            if x(1)>par_min
                x = [par_min,x];
                y = [min_val+threshold,y];
            end
            if x(end)<par_max
                x = [x,par_max];
                y = [y,min_val+threshold];
            end
        end
        plot(ax(i,j),x,y);
        yline(ax(i,j),min_val+threshold,"LineStyle","--")
        ax(i,j).YLim(2) = max(ax(i,j).YLim(2),min_val+1.1*threshold);
    end
end

for j = 1:npars
    title(ax(1,j),sm_par_display_names(j))
end
