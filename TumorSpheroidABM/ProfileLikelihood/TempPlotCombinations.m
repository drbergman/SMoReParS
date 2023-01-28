% a quick script for plotting parameter combinations after doing profile
% likelihood of the surrogate model on exp data (main_data.m)

figure;
ci=0;
for i = 1:2
    for j = i+1:3
        ci = ci+1;
        subplot(2,3,r2c(2,3,[1,ci]))
        plot(out{i}(i,:),out{i}(j,:))
        xlabel(para_names{i});ylabel(para_names{j})
        subplot(2,3,r2c(2,3,[2,ci]))
        plot(out{j}(j,:),out{j}(i,:))
        xlabel(para_names{j});ylabel(para_names{i})
    end
end