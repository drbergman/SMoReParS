% This script will compare fitting the ABM using SMoRe ParS against using
% the data directly.

clearvars;
addpath("~/Documents/MATLAB/myfunctions/")

save_fig_opts.save_figs = true;
save_fig_opts.reprint = false;
save_fig_opts.file_types = ["fig","png"];

name_suffix = "LMS";
admission_method = "all_profiles_resampled";
n_conditions = 3;
n_time_series = 2;
experimental_data_file = "../ODEFitting/data/ExperimentalData.mat";

load("data/AdmittedParameters_" + name_suffix + "_" + admission_method + ".mat","admitted_parameters","cohort_name","files")

% %% load cohort data and experimental data
C = load(sprintf("../../data/%s/output.mat",cohort_name),"ids","cohort_size","lattice_parameters","nsamps_per_condition");
A = load(files.abm_data_file,"D");
A.D = reshape(A.D,n_conditions,[]);
E = load(experimental_data_file,"D","t","C");
E.t(1) = [];
for ci = 1:n_conditions
    E.D(ci).A(1,:) = [];
    E.D(ci).S(1,:) = [];
end

nt = length(E.t);

%% Compute RSS values
RSS = zeros([n_time_series,n_conditions,numel(admitted_parameters)]);

for i = 1:numel(admitted_parameters)
    sim_data = arrayify(A.D(:,i),"A",1);
    sim_data = permute(sim_data,[1,3,2]);
    RSS(:,:,i) = computeRSSBreakdown(E,sim_data);
end

experimental_sd = arrayify(E.D,"S");
LL = -0.5*(nt*log(2*pi()) + sum(experimental_sd.^2,"all","omitnan") + squeeze(sum(RSS,1:2,"omitnan")));

%% stratify LL by SMoRe ParS Admission/Rejection
LL_admitted = LL(admitted_parameters(:));
LL_rejected = LL(~admitted_parameters(:));

%% plot histograms of both
f=figure("Name","LogLikeComparison_LMS"); hold on;
ax = gca;
[~,binEdges] = histcounts(LL);
histogram(LL_admitted,"BinEdges",binEdges,"FaceColor","green","Normalization","pdf","DisplayName","Selected Samples","EdgeColor","none")
histogram(LL_rejected,"BinEdges",binEdges,"FaceColor","red","Normalization","pdf","DisplayName","Rejected Samples","EdgeColor","none")
legend(ax,"Location","northwest","FontSize",16)
xlabel("log-likelihood","FontSize",16)
ylabel("PDF","FontSize",16)
saveFigures(f,save_fig_opts)

%% distribution relative likelihood to best setup
N = sum(admitted_parameters,"all");
divider_log = log(1e-5);
minRelProb_log = min(LL)-max(LL);
lbase_log = divider_log-minRelProb_log;
fn = @logLinScaleX;
LL_sort = sort(LL,"ascend");
XData = {LL_sort-max(LL),sort(LL_admitted)-max(LL),sort(LL_rejected)-max(LL)};
is_log_val = true;
CDF = {linspace(0,1,numel(LL)),linspace(0,1,length(LL_admitted)),linspace(0,1,length(LL_rejected))};
Names = ["All","Admitted","Rejected"];
Colors = ["b","g","r"];

%% plot cdf
f=figure("Name","RelativeLikelihoodCDF_LMS");
hold on;
for i = 1:3
    plot(fn(XData{i},divider_log,lbase_log,is_log_val),CDF{i},"DisplayName",Names(i),"Color",Colors(i),"LineWidth",2)
end
plot(fn(LL_sort(1:end-N)-LL_sort(end),divider_log,lbase_log,is_log_val),linspace(0,1,length(LL_rejected)),"DisplayName","Least Likely","Color",0.5*[1 1 1],"LineWidth",2,"LineStyle",":")
plot(fn(LL_sort(end-N+1:end)-LL_sort(end),divider_log,lbase_log,is_log_val),linspace(0,1,length(LL_admitted)),"DisplayName","Most Likely","Color",0*[1 1 1],"LineWidth",2,"LineStyle",":")
xlabel("Relative Likelihood","FontSize",16)
ylabel("CDF","FontSize",16)
legend(gca,"location","west","FontSize",16,"AutoUpdate","off")
xL = xlim;
xL(1) = fn(minRelProb_log,divider_log,lbase_log,is_log_val);
minPow10 = minRelProb_log*log10(exp(1));
maxPow10 = divider_log*log10(exp(1));
lowerPow10 = flip(maxPow10:ceil(0.5*(minPow10-maxPow10)):minPow10);
lowerXTicks = fn(lowerPow10.*log(10),divider_log,lbase_log,is_log_val);
upperXTicks = .2:.2:1;
xticks(unique([lowerXTicks,upperXTicks]))
lowerXTickLabels = strings(1,numel(lowerPow10));
for i = 1:numel(lowerPow10)
    temp_pow10 = round(lowerPow10(i));
    if 10^temp_pow10==0
        lowerXTickLabels(i) = sprintf("10^{%d}",temp_pow10);
    else
        lowerXTickLabels(i) = string(10^temp_pow10);
    end
end
xticklabels([lowerXTickLabels,string(upperXTicks)])
xL(2) = 1;
xlim(xL)
xline(exp(divider_log),"LineStyle","--")
annotation("textarrow",[.2875,.1875],[0.9,0.9],"String","Log Scale");
annotation("arrow",[0.3675,0.4675],[0.9,0.9]);
annotation("textarrow",[.7625,.8625],[0.15,0.15],"String","Linear Scale");
annotation("arrow",[0.6675,0.5675],[0.15,0.15]);

saveFigures(f,save_fig_opts)

%% attempt to plot PDF of above
f = figure;
hold on;
for i = 1:3
    y = diff(CDF{i}(:))./diff(exp(XData{i}(:)));
    plot(fn(.5*(XData{i}(1:end-1)+XData{i}(2:end)),divider_log,lbase_log,is_log_val),y/max(y),"DisplayName",Names(i),"Color",Colors(i),"LineWidth",2)
end
set(gca,"YScale","log")
xlabel("Relative Likelihood","FontSize",16)
ylabel("PDF","FontSize",16)
legend(gca,"location","best","FontSize",16,"AutoUpdate","off")

%% attempt #2 to plot PDF of above
f = figure("Name","RelativeLikelihoodPDF_LMS");
hold on;
% x = [0,logspace(-30,0,20)];
x = linspace(0,1,26);
for i = 1:3
    if is_log_val
        temp = histcounts(XData{i},log(x),"Normalization","pdf");
        plot(.5*(x(1:end-1)+x(2:end)),temp,"DisplayName",Names(i),"Color",Colors(i),"LineWidth",2)
    else
        temp = histcounts(XData{i},x,"Normalization","pdf");
        plot(.5*(x(1:end-1)+x(2:end)),temp,"DisplayName",Names(i),"Color",Colors(i),"LineWidth",2)
    end
end
set(gca,"YScale","log")
xlabel("Relative Likelihood","FontSize",16)
ylabel("PDF","FontSize",16)
legend(gca,"location","best","FontSize",16,"AutoUpdate","off")
saveFigures(f,save_fig_opts)

%% plotting percentage of N best LLs that were admitted
f=figure("Name","AdmissionOfMostLikely_LMS");
[~,LL_order_descend] = sort(LL,"descend");
admitted_sort = admitted_parameters(LL_order_descend);
y = zeros(2);
y(1,1) = sum(admitted_sort(1:N))/N; % percentage of admitted that are in the most likely category
y(2,1) = sum(~admitted_sort(1:N))/(numel(LL)-N); % percentage of rejected that are in the most likely category
y(1,2) = sum(admitted_sort(N+1:end))/N; % percentage of admitted that are in the most likely category
y(2,2) = sum(~admitted_sort(N+1:end))/(numel(LL)-N); % percentage of rejected that are in the most likely category
% y(:,1) = [sum(admitted_sort(1:N)),sum(~admitted_sort(1:N))]/N;
% y(:,2) = [sum(admitted_sort(N+1:end)),sum(~admitted_sort(N+1:end))]/(numel(LL)-N);
bar(categorical(["Admitted","Rejected"],["Admitted","Rejected"]),100*y,"stacked")
ylabel("Percentage")
set(gca,"FontSize",16)
title("SMoRe ParS Admission of Most Likely Parameters")
legend(["Most Likely","Least Likely"],"Location","northwest")
saveFigures(f,save_fig_opts)

%% reorient above to have categories be most likely and least likely and compare 
f = figure("Name","AdmissionByLikelihoodQuantile_LMS");
ax=gca;
[LL_sort,LL_order_ascend] = sort(LL,"ascend");
admitted_sort = admitted_parameters(LL_order_ascend);
n_quants = 20;
% Q = quantile(LL_sort,[1-N/numel(LL),1]);
Q = quantile(LL_sort,n_quants);
I_prev = 1;
y = zeros(n_quants,2);
for i = 1:n_quants
    I = find(LL_sort<=Q(i),1,"last");
    y(i,1) = sum(admitted_sort(I_prev:I))/(I-I_prev+1);
    y(i,2) = 1-y(i,1);
    I_prev = I;
end
bar(linspace(0,100,n_quants),100*y,"stacked","BarWidth",1,"EdgeColor","none");
xlim([0,100]+50/(n_quants-1)*[-1 1])
ylim([0 100])
legend(["Admitted","Rejected"],"location","northwest","AutoUpdate","off")
ylabel("Probability")
xlabel("Likelihood Percentile")
set(ax,"FontSize",16)
yticks(ax,0:25:100)
ax.XAxis.TickLength = [0 0];
xline(100*(1-N/numel(LL)))

saveFigures(f,save_fig_opts)

%% Residuals of best pars by LL
[~,LL_order_descend] = sort(LL,"descend");

for i = 1:numel(C.lattice_parameters)
    if iscell(C.lattice_parameters(i).path)
        continue
    end
    if all(C.lattice_parameters(i).path==["chemo_pars","concentration"])
        condition_dim = i;
        break;
    end
end
parameter_dims = setdiff(1:length(C.lattice_parameters),condition_dim);
sample_dim = length(C.lattice_parameters)+1;

C.ids = permute(C.ids,[condition_dim,parameter_dims,sample_dim]);
C.ids = reshape(C.ids,n_conditions,[]);

%% compute all Z scores
Z_scores_all = zeros([nt,n_time_series,size(C.ids)]);
for i = 1:size(C.ids,2)
    for ci = 1:n_conditions
        sim_data = zeros(nt,n_time_series);
        load(sprintf("../../data/sims/%s/output_final.mat",C.ids(ci,i)),"tracked")
        if isfield(tracked,"NT")
            y_temp = tracked.NT;
        else
            y_temp = sum(tracked.phases,2);
        end
        total = sum(tracked.phases,2);
        state2_prop = sum(tracked.phases(:,[3,4,7,8]),2)./total;
        z_score = ([total(2:end),state2_prop(2:end)]-E.D(ci).A)./E.D(ci).S;
        if any(isnan(z_score),"all") && (ci~=1 || any(isnan(z_score(:,1))))
            disp('')
        end
        Z_scores_all(:,:,ci,i) = z_score;
    end
end
Z_scores_all = reshape(Z_scores_all,nt,n_time_series,n_conditions,length(LL_order_descend),C.nsamps_per_condition);

%% Z scores by likelihood

Z_scores = {zeros(nt,n_time_series,n_conditions,N,C.nsamps_per_condition),zeros(nt,n_time_series,n_conditions,numel(LL)-N,C.nsamps_per_condition)};
for i = 1:length(LL_order_descend)
    if i<=N
        Z_scores{1}(:,:,:,i,:) = Z_scores_all(:,:,:,LL_order_descend(i),:);
    else
        Z_scores{2}(:,:,:,i-N,:) = Z_scores_all(:,:,:,LL_order_descend(i),:);
    end
end

%% make the figure of these residuals
f = figureOnRight("Name","TimeSeriesZScores_ByLikelihood_LMS");
nr = n_time_series*n_conditions - 1; % number of rows (do not do a row for G2/M in control)
ax = gobjects(nr,nt);
for ti = 1:nt
    for ri = 1:nr
        ind = ri + (ri>=2); % skip what would be the second, i.e. G2/M proportion in control
        [tsi,ci] = ind2sub([n_time_series,n_conditions],ind);
        if isnan(E.D(ci).A(ti,tsi)) || isnan(E.D(ci).S(ti,tsi)) % if the data (or its SD) is NaN, then skip
            continue
        end
        min_temp = min(min(Z_scores{1}(ti,tsi,ci,:)),min(Z_scores{2}(ti,tsi,ci,:))); % minimum z_score at time tt(ti)
        max_temp = max(max(Z_scores{1}(ti,tsi,ci,:)),max(Z_scores{2}(ti,tsi,ci,:))); % maximum z_score at time tt(ti)

        ax(ri,ti) = subplot(nr,nt,r2c(nr,nt,[ri,ti])); hold on;
        [~,binEdges] = histcounts(cat(4,Z_scores{1}(ti,tsi,ci,:),Z_scores{2}(ti,tsi,ci,:)));
        histogram(ax(ri,ti),Z_scores{1}(ti,tsi,ci,:),"BinEdges",binEdges,"FaceColor","green","Normalization","pdf","DisplayName","Most Likely","EdgeColor","none")
        histogram(ax(ri,ti),Z_scores{2}(ti,tsi,ci,:),"BinEdges",binEdges,"FaceColor","red","Normalization","pdf","DisplayName","Least Likely","EdgeColor","none")
        if min_temp~=0 || max_temp~=0
            xlim(max(abs([min_temp,max_temp])) * [-1 1]) % set xlim to be big enough to get largest |residual| and center at x=0
        end
        if ri==1
            title(ax(ri,ti),sprintf("t = %3.2f d",E.t(ti)))
        end
        if ri==nr
            xlabel(ax(ri,ti),"Z Score")
        end
        if ti==1
            switch ri
                case 1
                    ylabel(["Control","Total"])
                case 2
                    ylabel(["Dose=0.75\muM","Total"])
                case 3
                    ylabel(["Dose=0.75\muM","G2/M Prop"])
                case 4
                    ylabel(["Dose=7.55\muM","Total"])
                case 5
                    ylabel(["Dose=7.55\muM","G2/M Prop"])
            end

        end
    end
end
legend(ax(1,1))
normalizeXLims(f)

saveFigures(f,save_fig_opts)

%% compare residuals by admittance-rejection

Z_scores_admitted = Z_scores_all(:,:,:,admitted_parameters,:);
Z_scores_rejected = Z_scores_all(:,:,:,~admitted_parameters,:);

%% make the figure of these residuals
f = figureOnRight("Name","TimeSeriesZScores_ByAdmittance_LMS");
nr = n_time_series*n_conditions - 1; % number of rows (do not do a row for G2/M in control)
ax = gobjects(nr,nt);
for ti = 1:nt
    for ri = 1:nr
        ind = ri + (ri>=2); % skip what would be the second, i.e. G2/M proportion in control
        [tsi,ci] = ind2sub([n_time_series,n_conditions],ind);
        if isnan(E.D(ci).A(ti,tsi)) || isnan(E.D(ci).S(ti,tsi)) % if the data (or its SD) is NaN, then skip
            continue
        end
        min_temp = min(min(Z_scores_admitted(ti,tsi,ci,:)),min(Z_scores_rejected(ti,tsi,ci,:))); % minimum z_score at time tt(ti)
        max_temp = max(max(Z_scores_admitted(ti,tsi,ci,:)),max(Z_scores_rejected(ti,tsi,ci,:))); % maximum z_score at time tt(ti)

        ax(ri,ti) = subplot(nr,nt,r2c(nr,nt,[ri,ti])); hold on;
        [~,binEdges] = histcounts(cat(4,Z_scores_admitted(ti,tsi,ci,:),Z_scores_rejected(ti,tsi,ci,:)));
        histogram(ax(ri,ti),Z_scores_admitted(ti,tsi,ci,:),"BinEdges",binEdges,"FaceColor","green","Normalization","pdf","DisplayName","Admitted","EdgeColor","none")
        histogram(ax(ri,ti),Z_scores_rejected(ti,tsi,ci,:),"BinEdges",binEdges,"FaceColor","red","Normalization","pdf","DisplayName","Rejected","EdgeColor","none")
        if min_temp~=0 || max_temp~=0
            xlim(max(abs([min_temp,max_temp])) * [-1 1]) % set xlim to be big enough to get largest |residual| and center at x=0
        end
        if ri==1
            title(ax(ri,ti),sprintf("t = %3.2f d",E.t(ti)))
        end
        if ri==nr
            xlabel(ax(ri,ti),"Z Score")
        end
        if ti==1
            switch ri
                case 1
                    ylabel(["Control","Total"])
                case 2
                    ylabel(["Dose=0.75\muM","Total"])
                case 3
                    ylabel(["Dose=0.75\muM","G2/M Prop"])
                case 4
                    ylabel(["Dose=7.55\muM","Total"])
                case 5
                    ylabel(["Dose=7.55\muM","G2/M Prop"])
            end

        end
    end
end
legend(ax(1,1))
normalizeXLims(f)

saveFigures(f,save_fig_opts)

%% set x data into log-lin scale
% x data is in log scale below the divider and linear scale above the divider
% the lbase is the logarithm used on the log scale side
function x = logLinScaleX(x,divider,lbase,is_log_val)

assert(isequal(x,sort(x))) % make sure x sorted data
if ~is_log_val
    assert(all(x>0)) % make sure all x>0 so the log plot makes sense
end

lowerX = x < divider;
if any(lowerX)
    maxLowerX = max(x(lowerX));
    if is_log_val
        x(lowerX) = x(lowerX)/lbase;
        c = exp(divider) - (divider-maxLowerX)/lbase;
    else
        x(lowerX) = log(x(lowerX))/log(lbase);
        c = divider - log(divider/maxLowerX)/log(lbase);
    end
    x(lowerX) = x(lowerX) + c - max(x(lowerX));
end

if is_log_val && any(~lowerX)
    x(~lowerX) = exp(x(~lowerX));
end

end