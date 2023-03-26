% fits the ODE parameters to the experimental data. data is scaled so the
% initial cell count is 100

clearvars;

p(1) = 24/11; % lambda
p(2) = 24/8;
p(3) = 24/4; % alpha
p(4) = 24/1;
p(5) = 1e3; % K
p(6) = 0.09; % chemo-induced death probability per uM of drug

lb = [0;0;0;0;0;0];
ub = [Inf;Inf;Inf;Inf;1e4;1/7.55];

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
        % F = @(p) sum(arrayfun(@(cci) sum(((computeTimeSeriesChemo(p,tt,doses(cci),chemo_death_is_continuous)-data(:,cci))./data_std(:,cci)).^2,'all'),1:3));

    % F = @(p) sum(arrayfun(@(cci) sum(((sum(computeTimeSeriesChemo(p,tt,doses(cci),chemo_death_is_continuous),2) - data(:,cci))./data_std(:,cci)).^2,'all'),1:3));

    % F = @(p) sum(((sum(computeTimeSeriesChemo(p,tt,doses(1),chemo_death_is_continuous),2) - data{1})./data_std{1}).^2,'all') + ...
    %     sum(arrayfun(@(cci) sum(((computeTimeSeriesChemo(p,tt,doses(cci),chemo_death_is_continuous) - data{cci})./data_std{cci}).^2,'all'),2:3));
    [pstar,fstar] = fmincon(@(p) F(p,tt,doses,data,data_std),x0,[],[],[],[],lb,ub,[],opts);

%%
figure;
tfull = linspace(0,3,100);
for i = 1:3
    ax(i) = subplot(1,3,i); hold on;
    yy = computeTimeSeries4Phases(pstar,tfull,doses(i));
    if i == 1
        fit_curve = plot(tfull,sum(yy,2),"--","LineWidth",2,"DisplayName","Fit");
    else
        plot(tfull,[sum(yy(:,1:2),2),sum(yy(:,3:4),2)],"--","LineWidth",2)
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

function out = F(p,tt,doses,data,data_std)

out = sum(((sum(computeTimeSeries4Phases(p,tt,doses(1)),2) - data{1})./data_std{1}).^2,'all');
for j = 2:3
    temp = computeTimeSeries4Phases(p,tt,doses(j));
    temp = [sum(temp(:,1:2),2),sum(temp(:,3:4),2)];
    out = out + sum(((temp-data{j})./data_std{j}).^2,'all');
end

end
