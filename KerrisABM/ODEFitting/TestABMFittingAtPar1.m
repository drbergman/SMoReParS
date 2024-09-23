clearvars

% 24-04-23: manually test the fitting of the three SM at ABM par #1


%% calculate AIC
model_type = ["exponential","logistic","von_bertalanffy"];

data_file = "../PostAnalysis/data/summary.mat";

n_models = numel(model_type);
load(data_file,"t","D","C","cohort_size")
AIC = zeros([n_models,prod(cohort_size)]);
for i = 1:n_models
    % if model_type(i)=="von_bertalanffy"
    %     model_ode_fitting = load(sprintf("data/OptimalParameters_%s_resampled.mat",model_type(i)),"fstar","resample_t");
    % else
    model_ode_fitting = load(sprintf("data/OptimalParameters_%s.mat",model_type(i)),"fstar","resample_t");
    % end
    % model_ode_fitting = load(sprintf("data/OptimalParameters_%s.mat",model_type(i)),"fstar","resample_t");
    if isfield(model_ode_fitting,"resample_t")
        model_ode_fitting.fstar = model_ode_fitting.fstar ./ length(model_ode_fitting.resample_t);
    else
        model_ode_fitting.fstar = model_ode_fitting.fstar ./ 300;
    end
    AIC(i,:) = 2*i + model_ode_fitting.fstar(:);
end

%%

ind_to_test = find(AIC(2,:) < AIC(3,:));

%%
for i = 1:length(ind_to_test)


    I = ind_to_test(i);
    exp_data = load("data/OptimalParameters_exponential.mat");
    log_data = load("data/OptimalParameters_logistic.mat");
    vB_data = load("data/OptimalParameters_von_bertalanffy.mat");

    fstars = [exp_data.fstar(I),log_data.fstar(I),vB_data.fstar(I)];

    fstars_normalized = fstars ./ [5,5,301];

    aic_this = fstars_normalized + 2*[1,2,3];

    rel_log = min(aic_this) - aic_this;
    probs = exp(0.5*rel_log);

    %% compute time series
    tt = 15:15:75;
    exp_sol = 100*exp(exp_data.P(I) * tt);

    p = log_data.P(:,I);
    log_sol = 100*p(2)./((p(2)-100).*exp(-p(1).*tt')+100);

    p = vB_data.P(:,I);
    fn_opts.model_type = "von_bertalanffy";
    vB_sol = computeTimeSeries(p,linspace(0,75,301),[],fn_opts,[]);

    %%
    files.data = "../PostAnalysis/data/summary.mat";
    D = load(files.data);
    data = D.D;
    data = data(I);

    %%
    figure; hold on;
    plot(tt,exp_sol)
    plot(tt,log_sol)
    plot(linspace(0,75,301),vB_sol)
    plot(linspace(0,75,301),data.A,"Color","black")
    plot(linspace(0,75,301),data.A-data.S,"Color","black","LineStyle","--")
    plot(linspace(0,75,301),data.A+data.S,"Color","black","LineStyle","--")

    %% rss
    tI = 61:60:301;
    rss_exp = sum(((exp_sol' - data.A(tI))./data.S(tI)).^2);
    rss_log = sum(((log_sol - data.A(tI))./data.S(tI)).^2);
    rss_vB = sum(((vB_sol(2:end) - data.A(2:end))./data.S(2:end)).^2);

end