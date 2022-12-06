clearvars;
addpath("~/Documents/MATLAB/myfunctions/")

cohort_name = "cohort_221012070831969";
load(sprintf("../data/%s/%s.mat",cohort_name,cohort_name))
fgfr3_effects = ["No effect","C (cytotoxic)";"R (recruit)","R+C"];
tumor_type = ["LA","HA";"LA + Mut","HA + Mut"];

if ~exist("nsamps_per_condition","var")
    nsamps_per_condition = numel(tracked);
end

%% make standard figures
options.print_everything = false;
if options.print_everything
    nfigs_per_event = 4 + 1 + 4 + 2;
else
    nfigs_per_event = 2;
end
f = gobjects(4*(nfigs_per_event),1);
f(1:nfigs_per_event) = plotEventByTumorType(tracked,"tum_apop",fgfr3_effects,nsamps_per_condition,"Apoptosis",options);
f((nfigs_per_event+1):2*nfigs_per_event) = plotEventByTumorType(tracked,"tum_prolif",fgfr3_effects,nsamps_per_condition,"Proliferation",options);
f((2*nfigs_per_event+1):3*nfigs_per_event) = plotEventByTumorType(tracked,"tum_contact_inhibition",fgfr3_effects,nsamps_per_condition,"Contact Inhibition",options);
f((3*nfigs_per_event+1):4*nfigs_per_event) = plotEventByTumorType(tracked,"imm_cleared",fgfr3_effects,nsamps_per_condition,"Immune Clearance",options);


%% make other figures
n_cohorts = numel(tracked)/nsamps_per_condition;
all_nr = ceil(sqrt(n_cohorts));
all_nc = ceil(n_cohorts/all_nr);
colors = lines(n_cohorts);
type_colors = jet(4);

g = gobjects(0,1);
t = arrayifyNonuniform(tracked,"t");
t = reshape(t,n_cohorts,nsamps_per_condition,[]);
TumProlif = arrayifyNonuniform(tracked,"tum_prolif");
TumProlif = reshape(TumProlif,n_cohorts,nsamps_per_condition,[],2,2);
TumContactInhibition = arrayifyNonuniform(tracked,"tum_contact_inhibition");
TumContactInhibition = reshape(TumContactInhibition,n_cohorts,nsamps_per_condition,[],2,2);
Weight = TumProlif + TumContactInhibition;
Proportion = TumProlif./Weight;
Avg = sum(Proportion.*Weight,2)./sum(Weight,2);
STD = zeros(size(Avg));
for ci = 1:n_cohorts
    for mut_ind = 1:2
        for ant_ind = 1:2
            for ti = 1:size(Proportion,3)
                w = Weight(ci,:,ti,mut_ind,ant_ind);
                w(isnan(w)) = 0;
                STD(ci,1,ti,mut_ind,ant_ind) = std(Proportion(ci,:,ti,mut_ind,ant_ind),w);
            end
        end
    end
end
max_t = max(t,[],'all');
h = 4/24;
n = ceil(max_t*4);

%% comparison of tumor proliferation vs contact inhibition rates to assess where in tumor these cell types are (more central should mean more contact inhibition)
g(end+1,1) = figure("Name","tum_prolif_contact_inhibition_comparison");
ax = gobjects(2,2);
pp = gobjects(n_cohorts,2,2);
ll = gobjects(n_cohorts,2,2);


for mut_ind = 1:2
    for ant_ind = 1:2
        ax(mut_ind,ant_ind) = subplot(2,2,r2c(2,2,[mut_ind,ant_ind]));
        hold on;
        for ci = 1:n_cohorts
            xx = squeeze(t(ci,1,:));
            yy = squeeze(Avg(ci,1,:,mut_ind,ant_ind));
            r_mean = ksr(xx,yy,h,n);
            pc{1} = [r_mean.x,flip(r_mean.x)];
            std_temp = squeeze(STD(ci,1,:,mut_ind,ant_ind));
            r_std = ksr(xx,std_temp,h,n);
            pc{2} = [r_mean.f - r_std.f,flip(r_mean.f+r_std.f)];
            pp(ci,mut_ind,ant_ind) = patch(ax(mut_ind,ant_ind),pc{1},min(1,max(0,pc{2})),colors(ci,:),'FaceAlpha',0.15,'EdgeColor','black');
            ll(ci,mut_ind,ant_ind) = plot(ax(mut_ind,ant_ind),r_mean.x,r_mean.f,'Color',colors(ci,:),"LineWidth",1.25);
        end
        title(ax(mut_ind,ant_ind),sprintf("%s %s",tumor_type(mut_ind,ant_ind),"Proliferation:Contact Inhibition as Probability"))
    end
end

set(ax,"XLim",[0,max_t])
normalizeXLims(g(end))
normalizeYLims(g(end))
g(end).Position(3:4) = [1263 516];
L = legend(ax(1,1),reshape(ll(:,1,1),[],1),fgfr3_effects(:),"Location","best");
set(ax,'FontSize',16)
xlabel(ax,'Time (days)')
ylabel(ax,sprintf("%s","Probability"))

%% comparison of tumor proliferation vs contact inhibition rates by condition to assess where in tumor these cell types are (more central should mean more contact inhibition)
g(end+1,1) = figure("Name","tum_prolif_contact_inhibition_comparison_by_condition");
ax = gobjects(n_cohorts,1);
pp = gobjects(n_cohorts,2,2);
ll = gobjects(n_cohorts,2,2);

for ci = 1:n_cohorts
    ax(ci) = subplot(all_nr,all_nc,r2c(all_nr,all_nc,ci));
    hold on;
    for mut_ind = 1:2
        for ant_ind = 1:2
            xx = squeeze(t(ci,1,:));
            yy = squeeze(Avg(ci,1,:,mut_ind,ant_ind));
            r_mean = ksr(xx,yy,h,n);
            pc{1} = [r_mean.x,flip(r_mean.x)];
            std_temp = squeeze(STD(ci,1,:,mut_ind,ant_ind));
            r_std = ksr(xx,std_temp,h,n);
            pc{2} = [r_mean.f - r_std.f,flip(r_mean.f+r_std.f)];
            pp(ci,mut_ind,ant_ind) = patch(ax(ci),pc{1},min(1,max(0,pc{2})),type_colors(mut_ind+(ant_ind-1)*2,:),'FaceAlpha',0.15,'EdgeColor','black');
            ll(ci,mut_ind,ant_ind) = plot(ax(ci),r_mean.x,r_mean.f,'Color',type_colors(mut_ind+(ant_ind-1)*2,:),"LineWidth",1.25);
        end
    end
    title(ax(ci),sprintf("%s %s",fgfr3_effects(ci),"Proliferation:Contact Inhibition as Probability"))
end

set(ax,"XLim",[0,max_t])
normalizeXLims(g(end))
normalizeYLims(g(end))
g(end).Position(3:4) = [1263 516];
L = legend(ax(1),reshape(ll(1,:,:),[],1),tumor_type(:),"Location","best");
set(ax,'FontSize',16)
xlabel(ax,'Time (days)')
ylabel(ax,sprintf("%s","Probability"))

%% print
printFigures(f,cohort_name);
printFigures(g,cohort_name);
