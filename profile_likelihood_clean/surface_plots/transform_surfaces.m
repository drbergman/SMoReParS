


fig = openfig("test_theta_abm_plot.fig");
plt = findobj(fig, 'type', 'surface')

X = plt.XData;
Y = plt.YData;
upper_bd = plt(1).ZData;
upper_bd = 1 - upper_bd.^(-1);
lower_bd = plt(2).ZData;
lower_bd = 1 - lower_bd.^(-1);
size(lower_bd)

figure(1)
surf(X, Y, lower_bd, 'FaceColor', 'b', 'FaceAlpha', 0.5, 'EdgeColor', 'none')
hold on;
surf(X,Y,upper_bd, 'FaceColor', 'm', 'FaceAlpha', 0.5, 'EdgeColor', 'none')
xlabel("ABM div_{lim}");
ylabel("ABM p_{div}");
zlabel("\gamma");
savefig("test_2_gamma_abm_plot.fig");


