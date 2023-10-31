clearvars;

save_fig_opts.save_figs = true;
save_fig_opts.reprint = true;
save_fig_opts.file_types = ["fig","png"];
save_fig_opts.resolution = '-r1200';

file_name = "ExperimentalData";

load("data/ExperimentalData.mat","D","t")

nr = 1;
nc = 2;

f = figureOnRight("Name",file_name);
ax = gobjects(nr,nc);
color = [113,191,68;234,231,35;238,49,36]/255;
for tsi = 1:2
    ax(tsi) = subplot(nr,nc,tsi);
    hold on;
    for ci = 1:3
        plot(ax(tsi),t,D(ci).A(:,tsi),"Color",color(ci,:),"LineWidth",0.5)
        errorbar(t,D(ci).A(:,tsi),D(ci).S(:,tsi),"vertical","Color",color(ci,:),"LineStyle","none","CapSize",2);
    end
end

xlabel(ax,"Time (d)")
ylabel(ax(1),"Total Cell Count")
ylabel(ax(2),"G2/M Fraction")

%%
f.Units = "inches";
f.Position(3) = 1.5;
f.Position(4) = .8;

set(ax,"FontSize",6)

%%
margin = struct("left",0.18,"right",.04,"top",0.09,"bottom",0.27);
spacing = struct("horizontal",0.22,"vertical",0.1);
uniformAxisSpacing(ax,margin,spacing);

%%
saveFigures(f,save_fig_opts)