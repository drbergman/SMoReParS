function out = moatSample(x,par_names,D,I,nsamps,alpha,ci_relative_spread)

% runs nsamps at the point in parameter space defined by x. M has the base
% parameters of the simulation. par_names stores the correspondence between
% indices in x and parameters in M. D is a vector of distributions for the
% parameters, so that a call to ICDF can return the appropriate value.
p = zeros(2,1);
for i = 1:numel(x)
    p(I(par_names(i))) = icdf(D(par_names(i)),x(i));
end

NT = complexModel(repmat(p,1,nsamps));

% std(NT)/sqrt(numel(NT)) gives the estimated SD of means of sample size numel(NT)
% tinv(1-0.5*alpha,numel(NT)-1) gives the t value for the 95% CI with numel(NT)-1 degrees of freedom
% Thus, the product of these (call it R) gives the radius of the CI about mean(NT) for estimating the true population mean, i.e. true ABM output on Day 3
% I want to insist that this radius, or spread, is less than 10% of the estimated mean, i.e. R/mean(NT) < 0.1, otherwise I continue taking more samples

relative_spread = tinv(1-0.5*alpha,numel(NT)-1)*std(NT)/(sqrt(numel(NT))*mean(NT));

if relative_spread > ci_relative_spread
    while relative_spread > ci_relative_spread % while the confidence interval for the mean is wider than 10% of the mean, keep going
        NT(end+1) = complexModel(p);
        relative_spread = tinv(1-0.5*alpha,numel(NT)-1)*std(NT)/(sqrt(numel(NT))*mean(NT));
%         fprintf("  Spread after %d: %3.2f%%\n",numel(NT),100*relative_spread)
    end
else
%     fprintf("Spread after 10: %3.2f%%\n",100*relative_spread)
end

out = mean(NT);

