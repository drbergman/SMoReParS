
clear all;
close all;

fig = openfig(".\no_confidence_line\theta_resnorm_plot.fig");
plt = findobj(fig, 'type', 'line')
X = plt.XData;
Y = plt.YData;
hold off;

%%{
X_new = 1 - X.^(-1);
figure(2)
plot(X_new, Y)
hold on;

X2 = 1 - [1.56313, 2.5].^(-1);
Y2 = [145.79, 145.79];
plot(X2, Y2, 'Color', 'r');

title("Residual Norm");
xlabel("\gamma");
ylabel("Residual Norm");
legend("Residual Norm", "Confidence Threshold");

savefig("gamma_resnorm_plot.fig");


%}


