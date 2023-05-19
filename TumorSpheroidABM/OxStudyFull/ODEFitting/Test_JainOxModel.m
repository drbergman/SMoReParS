% a script to test the Jain Oxaliplatin model

clearvars;

save_fig_opts.file_types = ["fig","png"];

lambda = 5.96;
alphaRP = 0.99;
theta = 1;
VT = 5.48; % in thousands of cells
V0 = .55; % in thousands of cells
alphaP = 8.69;
kalpha = 21.67;
a = 0.5;
rho0 = 0.06; % not listed in their Table 1 under Oxaliplatin, using the value reported in Taxol model of Table 1
delta0 = 0.96;
kdelta = 9.22;
b = 1;
alphaR = 3.0;
p = [lambda;alphaRP;theta;VT;V0;alphaP;kalpha;a;rho0;delta0;kdelta;b;alphaR];

IC = [.10;.90;0;0];
dose = [0;0.75;7.55];

tt = linspace(0,3,100);

D = load("data/ExperimentalData.mat","D","t");

%% plotting
f=figure("Name","JainOxModelFit","Position",[0 0 1000 420]);
ax = gobjects(2,1);
for j = 1:2
    ax(j) = subplot(1,2,j);
    hold on;
end
title(ax(1),"Total Cell Count")
title(ax(2),"G2/M Proportion")
colors = lines(3);
for i = 1:3
    out = computeTimeSeries_jain(p,tt,dose(i));

    ll(i) = plot(ax(1),tt,out(:,1),"LineWidth",2,"DisplayName",num2str(dose(i)) + "\muM","Color",colors(i,:));
    plot(ax(2),tt,out(:,2),"LineWidth",2,"DisplayName",num2str(dose(i)) + "\muM","Color",colors(i,:))

    scatter(ax(1),D.t,D.D(i).A(:,1)/100,"filled","MarkerFaceColor",colors(i,:),"Marker","s","DisplayName","Data")
    errorbar(ax(1),D.t,D.D(i).A(:,1)/100,D.D(i).S(:,1)/100,"Color",colors(i,:),"LineStyle","none","LineWidth",1)
    scatter(ax(2),D.t,D.D(i).A(:,2),"filled","MarkerFaceColor",colors(i,:),"Marker","s","DisplayName","Data")
    errorbar(ax(2),D.t,D.D(i).A(:,2),D.D(i).S(:,2),"Color",colors(i,:),"LineStyle","none","LineWidth",1)

end
ylabel(ax(1),"Count in thousands")
ylabel(ax(2),"Proportion")
% ylim(ax(2),[0 1])
L=legend(ax(1),ll,"location","best");
L.Title.String = "Dose";
xlabel(ax,"Time (d)")
set(ax,"FontSize",16)

saveFigures(f,save_fig_opts)
