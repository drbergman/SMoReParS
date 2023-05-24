% A quick script to plot the state variables for one of the ABM parameter
% vectors

clearvars;
cohort_name = "cohort_230124175743017";
addpath("~/Documents/MATLAB/myfunctions/")
save_figs = false;

%% load data
load(sprintf("../../data/%s/output.mat",cohort_name),"ids","nsamps_per_condition","cohort_size","lattice_parameters");
nsamps_per_parameter_vector = nsamps_per_condition;

ids = reshape(ids,[],nsamps_per_parameter_vector);

I = randi(size(ids,1));

for i = nsamps_per_parameter_vector:-1:1
    S = load(sprintf("../../data/sims/%s/output_final.mat",ids(I,i)));
    if size(S.tracked.phases,2)>4
        error("Not sure how we will count the arrested compartment in this.")
    end
    phase_count(:,:,:,i) = reshape(S.tracked.phases,[],2,2);
end


t_abm = S.tracked.t;
t_abm = round(1440*t_abm)/1440; % to make sure that the last time point is actually 3 days (not 3-eps() days)
nt_abm = length(t_abm);

t = [0;10;24;48;72]/24;
nt = length(t);

phase_count = reshape(phase_count,[nt_abm,2,2,nsamps_per_parameter_vector]);
phase_count_sampled = interp1(t_abm,phase_count,t);
state_vars = sum(phase_count_sampled,2);
state_vars = reshape(state_vars,nt,2,[],nsamps_per_parameter_vector);
state_vars = permute(state_vars,[4,2,1,3]);

%% plot joint probabilities
nr = ceil(sqrt(nt));
nc = ceil(nt/nr);
f=figure;
th = linspace(0,2*pi(),1001);
st = sin(th);
ct = cos(th);
ax = gobjects(nt,1);
sds = sqrt(chi2inv(.95,2));
colors = lines(2);
for ti = 1:nt
    ax(ti) = subplot(nr,nc,ti);
    hold on
    X = state_vars(:,:,ti,1);
    C = cov(X);
    mu = mean(X);
    [eV,ev] = eig(C);
    d = sqrt(max(0,diag(ev))); % don't let it be negative because that only happens when a rounding error puts it just below 0
    % patch(xx(1,:)+mu(1),xx(2,:)+mu(2),[0.4,0.4,0.4])
    for j = 1:numel(sds)
        xx = eV*(sds(j)*[ct;st].*d);
        JJ(j) = plot(xx(1,:)+mu(1),xx(2,:)+mu(2),"Color",colors(1,:),"LineWidth",1,"DisplayName","Joint Probability");
    end
    for j = 1:numel(sds)
        yy = sds(j)*[ct;st].*sqrt(max(0,diag(C)));
        II(j) = plot(yy(1,:)+mu(1),yy(2,:)+mu(2),"Color",colors(2,:),"LineWidth",1,"DisplayName","Independent");
    end
    S = scatter(X(:,1),X(:,2),40,'k','filled',"DisplayName","Data");
    title(sprintf("t = %dh",t(ti)*24))
    % scatter(mu(1),mu(2),'filled')
    % axis square
end
% xlabel(ax,"G1/S")
% ylabel(ax,"G2/M")
set(ax,"FontSize",16)
% axis(ax,"square")
L = legend(ax(end),[JJ(1),II(1),S],"Location","best");
L.Position = [0.5830    0.1228    0.2598    0.1369];
%% save plots
if save_figs
    fig_names = "JointProbabilities";
    for i = 1:numel(f)
        if isempty(f(i).Name)
            f(i).Name = fig_names(i);
        end
        fig_folders = ["figures/fig","figures/png"];
        for j = 1:numel(fig_folders)
            if ~exist(fig_folders(j),"dir")
                mkdir(fig_folders(j))
            end
        end
        savefig(f(i),sprintf("figures/fig/%s",f(i).Name))
        print(f(i),sprintf("figures/png/%s",f(i).Name),"-dpng")
    end
end

