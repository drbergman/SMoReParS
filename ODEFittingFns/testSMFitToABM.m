function f = testSMFitToABM(par_file,data_file,nsamps,fn,fn_opts,par_names)

if contains(path,"myfunctions")
    path_changed = false;
else
    path_changed = true;
    addpath("~/Documents/MATLAB/myfunctions/")
end

load(par_file,"P")
load(data_file,"t","D","C","cohort_size","nsamps_per_parameter_vector","n_time_series","n_conditions"); % some of these variables are not used now, but they might be once I get to filling out the conditional statements below

P = reshape(P,size(P,1),[]);
I = randperm(size(P,2),nsamps);

f = gobjects(2,1);
%% these will be filled out as new cases are covered with this
if n_time_series == 1 && n_conditions == 1 % just plot in a rough square
    f(1)=figure;
    nr = ceil(sqrt(nsamps));
    nc = ceil(nsamps/nr);
    ax = gobjects(nsamps,1);
    for i = 1:nsamps
        ax(i) = subplot(nr,nc,i,"NextPlot","add");
        patch(ax(i),[t(:);flip(t(:))],[D(I(i)).A-D(I(i)).S;flipud(D(I(i)).A+D(I(i)).S)],"black","FaceAlpha",0.2,"EdgeColor","none")
        plot(ax(i),t,D(I(i)).A,"black")
        out = fn(P(:,I(i)),t,C{1},fn_opts);
        plot(ax(i),t,out,"--","LineWidth",2)
    end
elseif n_time_series > 1 && n_conditions == 1
    f(1)=figure;
    nr = nsamps;
    nc = n_time_series;
    ax = gobjects(nsamps,n_time_series);
    for ri = 1:nsamps
        out = fn(P(:,I(ri)),t,C{1},fn_opts);
        for ci = 1:n_time_series
            ax(ri,ci) = subplot(nr,nc,r2c(nr,nc,[ri,ci]),"NextPlot","add");
            patch(ax(ri,ci),[t(:);flip(t(:))],[D(I(ri)).A(:,ci)-D(I(ri)).S(:,ci);flipud(D(I(ri)).A(:,ci)+D(I(ri)).S(:,ci))],"black","FaceAlpha",0.2,"EdgeColor","none")
            plot(ax(ri,ci),t,D(I(ri)).A(:,ci),"black")
            plot(ax(ri,ci),t,out(:,ci),"--","LineWidth",2)
        end
    end
end
xlim(ax,[t(1) t(end)])
%% histograms of parameter values
f(2)=figure;
for i = 1:size(P,1)
    subplot(size(P,1),1,i)
    histogram(P(i,:))
    title(par_names(i))
end

if path_changed
    rmpath("~/Documents/MATLAB/myfunctions/")
end
