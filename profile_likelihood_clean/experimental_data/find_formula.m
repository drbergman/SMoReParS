

fig = openfig("mu_p_vs_theta.fig");
plt = findobj(fig, 'type', 'scatter')
mu = plt.XData;
mu = mu(1:81);
theta = plt.YData;
theta = theta(1:81);

f = @(p, xdata) p(1)*xdata + p(2);
[pfit, resnorm] = lsqcurvefit(f, [0, 0], mu, theta);

figure(2)
plot(mu, theta)
hold on
plot(mu, pfit(1)*mu + pfit(2))