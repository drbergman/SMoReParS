clearvars

save_fig_opts.save_figs = true;
save_fig_opts.reprint = false;
save_fig_opts.file_types = ["fig","png"];
save_fig_opts.resolution = '-r1200';

% fac = .8;
% colors = {"r","b","g"};
colors = {[237,35,39]/255,[108,190,70]/255,[58,84,164]/255};
% colors = {fac*[0,1,1],fac*[1,0,1],fac*[1,1,0]};
par_color = containers.Map(["\lambda","\alpha","K"],colors);

f = openfig("../figures/fig/ParTriplesAll_SMFitToData_New.fig");

for i = 1:numel(f.Children)
    if isa(f.Children(i),'matlab.graphics.illustration.Legend')
        delete(f.Children(i))
        break
    end
end

ax = gca;
for i = 1:numel(ax.Children)
    ax.Children(i).LineWidth = 0.5;
    ax.Children(i).LineStyle = "-";
    ax.Children(i).Marker = ".";
    ax.Children(i).MarkerSize = 4;
    ax.Children(i).Color = par_color(ax.Children(i).DisplayName);
end

ax.XAxis.Label.VerticalAlignment = "middle";
ax.YAxis.Label.VerticalAlignment = "middle";
ax.ZAxis.Label.VerticalAlignment = "middle";

%%

f.Units = "inches";
f.Position(3) = 2;
f.Position(4) = 2;

margin = struct("left",[],"right",[],"top",0,"bottom",.15);
spacing = struct("horizontal",0,"vertical",0);
uniformAxisSpacing(gca,margin,spacing);
view(ax,-33,19);
cd("..")
saveFigures(f,save_fig_opts)
cd("FigureScripts")

