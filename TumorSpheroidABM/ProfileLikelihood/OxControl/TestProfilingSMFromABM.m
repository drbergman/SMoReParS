clearvars;

load("data/ProfileLikelihoods.mat","out")

nsamps = 5;
I = randsample(2187,nsamps,false);
nr = size(out,1);
threshold = chi2inv(0.95,3);
figure;
for i = 1:nsamps
    for j = 1:nr
        subplot(nr,nsamps,r2c(nr,nsamps,[j,i]))
        plot(out{j,I(i)}(1,:),out{j,I(i)}(2,:))
        min_val = min(out{j,I(i)}(2,:));
        yline(min_val+threshold,"LineStyle","--")
    end
end