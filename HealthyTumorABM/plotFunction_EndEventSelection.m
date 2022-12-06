function plotFunction_EndEventSelection(M, prob_matrix)

%% healthy probabilities
order = [2,1]; % first, separate by fast vs slow, then engaged/non-engaged
xs = linspace(0,1,size(prob_matrix.healthy,1));
ys = sortrows(prob_matrix.healthy(:,order));
for i = 1:2
    set(M.fig.healthy_prob_bar(i),'XData',xs,'YData',ys(:,i))
end
if length(xs)>1
    M.fig.ax(M.fig.healthy_prob_ind).XLim = [0,1] + .5*[-1,1]*xs(2);
end

%% tumor probabilities
order = [2,1]; % first, separate by fast vs slow, then engaged/non-engaged
xs = linspace(0,1,size(prob_matrix.tumor,1));
ys = sortrows(prob_matrix.tumor(:,order));
for i = 1:2
    set(M.fig.tumor_prob_bar(i),'XData',xs,'YData',ys(:,i))
end
if length(xs)>1
    M.fig.ax(M.fig.tumor_prob_ind).XLim = [0,1] + .5*[-1,1]*xs(2);
end

