% a script to visualize (projections of) the hypersurfaces

clearvars;

cohort_name = "cohort_230124175743017";

load("data/ProfileLikelihoods.mat")
C = load(sprintf("../../data/%s/output.mat",cohort_name),"cohort_size","lattice_parameters");
threshold = chi2inv(0.95,3);
%%
out = reshape(out,[size(out,1),C.cohort_size]);

%% make the meshes
temp = out(:,:,:,1); % take a small slice of out to focus on for drawing sample surfaces

xx = C.lattice_parameters(1).values;
yy = C.lattice_parameters(2).values;

S_min = zeros(3,3,3);
S_max = zeros(3,3,3);
for i = 1:3 % an ABM parameter index
    for j = 1:3 % another ABM parameter index
        for pi = 1:3 % SM parameter index

            [S_min(i,j,pi),S_max(i,j,pi)] = getProfileBounds(temp{pi,i,j},threshold);


            % S_min(i,j,pi) = min(temp{pi,i,j}(1,:));
            % S_max(i,j,pi) = max(temp{pi,i,j}(1,:));
        end
    end
end

%% surface for all three ODE parameters projected onto first two abm parameter dimensions
figure;
for pi = 1:3
    ax(pi) = subplot(3,1,pi); hold on
    mesh(xx,yy,S_min(:,:,pi)',"FaceColor","blue","FaceAlpha",0.2,"EdgeColor","none")
    mesh(xx,yy,S_max(:,:,pi)',"FaceColor","red","FaceAlpha",0.2,"EdgeColor","none")
    view(3)
    xlabel("ABM Par 1")
    ylabel("ABM Par 2")
end
%% surface for K projected onto first two abm parameter dimensions
ode_par_name = {'\lambda','\alpha','K'};
file_name = {'lambda','alpha','K'};
for pi = 1:3
    figure;
    hold on
    mesh(xx,yy,S_min(:,:,pi)',"FaceColor","blue","FaceAlpha",0.2,"EdgeColor","none")
    mesh(xx,yy,S_max(:,:,pi)',"FaceColor","red","FaceAlpha",0.2,"EdgeColor","none")
    view(3)
    xlabel("ABM Par 1")
    ylabel("ABM Par 2")
    title(ode_par_name{pi})
    set(gca,'FontSize',20)
    savefig(sprintf("figs/sample_surface_%s.fig",file_name{pi}))
    print(sprintf("figs/sample_surface_%s.png",file_name{pi}),"-dpng")
end


%% a more refined mesh for the three surface (projected onto the first two abm parameter dimensions)
figure;
xxq = linspace(xx(1),xx(end),15);
yyq = linspace(yy(1),yy(end),15);
pars = {C.lattice_parameters.values};
pars = pars(1:2);
for pi = 1:3
    for i = 1:numel(xxq)
        for j = 1:numel(yyq)
            V_min(i,j,pi) = interpn(pars{:},S_min(:,:,pi),xxq(i),yyq(j),"linear");
            V_max(i,j,pi) = interpn(pars{:},S_max(:,:,pi),xxq(i),yyq(j),"linear");
        end
    end
    ax(pi) = subplot(3,1,pi); hold on;
    mesh(ax(pi),xxq,yyq,V_min(:,:,pi)',"FaceColor","blue","FaceAlpha",0.2,"EdgeColor","none")
    mesh(ax(pi),xxq,yyq,V_max(:,:,pi)',"FaceColor","red","FaceAlpha",0.2,"EdgeColor","none")
    view(3)
end
