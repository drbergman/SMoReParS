clearvars;

cohort_name = "cohort_230124175743017";

load("OptimalParameters_noapop.mat")
Sum = load(sprintf("../data/%s/summary.mat",cohort_name),"ode_state*");
load(sprintf("../data/%s/output.mat",cohort_name),"ids");
load(sprintf("../data/sims/%s/output_final.mat",ids(1)),"tracked");
tt = tracked.t;
I = randperm(numel(P)/size(P,1),5);

sz = size(ids);
sz(end) = [];
Sum.ode_state_average = reshape(Sum.ode_state_average,length(tt),2,[]);
Sum.ode_state_std = reshape(Sum.ode_state_std,length(tt),2,[]);
P = reshape(P,size(P,1),[]);

figure;
for i = 1:5
    subplot(5,1,i)
    hold on

    avg = Sum.ode_state_average(:,:,I(i));
    std = Sum.ode_state_std(:,:,I(i));
    for j = 1:2
        patch([tt;flip(tt)],[avg(:,j)-std(:,j);flipud(avg(:,j)+std(:,j))],"black","FaceAlpha",0.2,"EdgeColor","none")
    end
    plot(tt,avg,"black")

    out = computeTimeSeries(P(:,I(i)),tt);
    plot(tt,out,"--","LineWidth",2)
end

%%
% par_names = ["lambda","alpha","K","delta","G1 prop_0"];
par_names = ["lambda","alpha","K"];
figure;
for i = 1:3
    subplot(3,1,i)
    histogram(P(i,:))
    title(par_names(i))
end