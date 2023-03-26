% fits the ODE parameters to the experimental data. data is scaled so the
% initial cell count is 100

clearvars;

p(1) = 24/19; % lambda
p(2) = 24/5; % alpha
p(3) = 1e3; % K
p(4) = 0.1; % chemo-induced death rate per uM of drug
chemo_death_is_continuous = false; % does chemo death occur over entirety of each phase (true)? Or is it a one-time event during a phase and so it happens at a higher rate during shorter phases (false)?

lb = [0;0;0;0];
ub = [Inf;Inf;1e4;2];

opts = optimset('Display','off','TolFun',1e-12,'TolX',1e-12);
%%
% fn = fieldnames(p);
% npars = numel(fn);
% x0 = zeros(npars,1);
% for i = 1:npars
%     x0(i) = p.(fn{i});
% end
x0 = p;
npars = length(p);


% Data from Jang et al. Cancer Res Treat 2002;34:372. Millions of cells.
% Data point for 5 hours taken out, since it is incommensurate.

tt = [0      10     24     36     48     72    ]';        % hours
tt = tt/24;                                             % days

% Control data
data{1} =     [0.899  1.340  1.633  2.408  3.557  5.583]';   % millions of cells
data_std{1} = [0.099  0.193  0.207  0.298  0.168  0.364]';   % millions of cells
% dose = 0.75 uM
prop_in_phase1 = [0.91;0.81;0.82;0.85;0.88;0.9];
data{2} =    [0.899    1     1.633  2.310  2.558  3.289]' .* (prop_in_phase1.*[1,-1] + [0,1]);   % millions of cells
data_std{2} = [0.099  0.193  0.207  0.440  0.850  0.440]' .* (prop_in_phase1.*[1,-1] + [0,1]);   % millions of cells
% dose = 7.55 uM
prop_in_phase1 = [0.91;0.81;0.82;0.80;0.79;0.79];
data{3} =    [0.899   0.920  1.150  1.000  1.000  0.960]' .* (prop_in_phase1.*[1,-1] + [0,1]);   % millions of cells
data_std{3} = [0.099  0.263  0.207  0.210  0.200  0.300]' .* (prop_in_phase1.*[1,-1] + [0,1]);   % millions of cells

doses = [0;0.75;7.55];
%% scale data
for i = 1:3
    factor = 100./sum(data{i}(1,:));
    data_std{i} = data_std{i} .* factor;
    data{i} = data{i} .* factor;
end
%%
P = zeros(npars,1);

%     F = @(p) sum(((computeTimeSeries(p,tt)./cell_count - 1)).^2,'all');
%     F = @(p) sum(((computeTimeSeries(p,tt) - data)).^2,'all');
        % F = @(p) sum(arrayfun(@(cci) sum(((computeTimeSeries(p,tt,doses(cci),chemo_death_is_continuous)-data(:,cci))./data_std(:,cci)).^2,'all'),1:3));

    % F = @(p) sum(arrayfun(@(cci) sum(((sum(computeTimeSeries(p,tt,doses(cci),chemo_death_is_continuous),2) - data(:,cci))./data_std(:,cci)).^2,'all'),1:3));

    F = @(p) sum(((sum(computeTimeSeries(p,tt,doses(1),chemo_death_is_continuous),2) - data{1})./data_std{1}).^2,'all') + ...
        sum(arrayfun(@(cci) sum(((computeTimeSeries(p,tt,doses(cci),chemo_death_is_continuous) - data{cci})./data_std{cci}).^2,'all'),2:3));
    [pstar,fstar] = fmincon(F,x0,[],[],[],[],lb,ub,[],opts);

%%
figure;
tfull = linspace(0,3,100);
for i = 1:3
    ax(i) = subplot(1,3,i); hold on;
    yy = computeTimeSeries(pstar,tfull,doses(i),chemo_death_is_continuous);
    if i == 1
        fit_curve = plot(tfull,sum(yy,2),"--","LineWidth",2,"DisplayName","Fit");
    else
        plot(tfull,yy,"--","LineWidth",2)
    end
    for j = 1:size(data{i},2)
        data_curve = plot(tt,data{i}(:,j),"black","Marker","o","MarkerFaceColor","black","DisplayName","Data");
        data_patch = patch([tt;flip(tt)],[data{i}(:,j)-data_std{i}(:,j);flip(data{i}(:,j)+data_std{i}(:,j))],"black","FaceAlpha",0.2,"EdgeColor","none","DisplayName","+/- SD");
    end
    title(sprintf("C = %3.2fuM",doses(i)),"Interpreter","none")
end
xlabel(ax,"Time (d)")
ylabel(ax,"Scaled Cell Count")

% legend([fit_curve;data_curve;data_patch],"Location","northwest","FontSize",22)

set(ax,"FontSize",20)
