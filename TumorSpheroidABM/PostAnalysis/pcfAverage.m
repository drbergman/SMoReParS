function [t,Avg] = pcfAverage(pcf,options)

Avg = zeros(length(options.rr),0);
weights = zeros(1,0);
t = [];
for si = 1:numel(pcf)
    t = unique([t,pcf(si).t]);
    szs = [size(Avg,2),size(pcf(si).avg,2)]; % number of time points for each so far
    if szs(1)<szs(2)
        weights = weights+1;
        Avg = ((si-1)*Avg + pcf(si).avg(:,1:szs(1)))./weights;
        weights(end+1:szs(2)) = 1;
        Avg(:,end+1:szs(2)) = pcf(si).avg(:,szs(1)+1:end);
    elseif szs(1)==szs(2)
        weights = weights+1;
        Avg = ((si-1)*Avg + pcf(si).avg)./weights;
    else
        weights(1:szs(2)) = weights(1:szs(2))+1;
        Avg(:,1:szs(2)) = ((si-1)*Avg(:,1:szs(2)) + pcf(si).avg)./weights(1:szs(2));
    end
end
