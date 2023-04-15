function [f,I] = testProfileSMFromABM(profile_file,nsamps,sm_par_display_names)

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
