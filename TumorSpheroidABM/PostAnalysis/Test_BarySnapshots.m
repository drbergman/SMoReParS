clearvars;

v_ha_mut = [0;0];
v_la_mut = [0.5;sqrt(3)/2];
v_la = [1;0];

v = [v_la,v_la_mut,v_ha_mut];

T = 50; % day 10



%%
cohort_name = "cohort_221021113614089";
path_to_cohort_file = sprintf("../data/%s/%s.mat",cohort_name,cohort_name);
load(path_to_cohort_file,"total_runs","ids","nsamps_per_condition")

A10 = zeros(total_runs,3);
for si = 1:total_runs
    sim_name = ids(si);
    path_to_sim_folder = sprintf("../data/%s",sim_name);
    load(sprintf("%s/output_final.mat",path_to_sim_folder),"tracked")
    temp = reshape(tracked.phase_cell_hours,[],4);
    temp = temp(:,[1,2,4]);
    if T<=tracked.t(end)
        temp = interp1(tracked.t,temp,T,"linear");
    else
        ind = find(any(temp>0,2),1,"last");
        temp = temp(ind,:);
    end
    temp = temp./sum(temp,2);
    A10(si,:) = temp;
end

A10 = reshape(A10,2,2,4,10,3);

%%

nbins = 11;
for fdi = 1:4 % fig dimensions index

    f(fdi)=figure;
    scatter(v(1,:),v(2,:),'filled')
    f(fdi).Position(1) = 1547;
    f(fdi).Position(3:4) = [664 543];
    for ri = 1:2
        for ci = 1:2
            ax(ri,ci) = subplot(2,2,r2c(2,2,[ri,ci]));
            [val,xe,ye] = histcounts2(A10(ri,ci,fdi,:,1),A10(ri,ci,fdi,:,2),linspace(0,1,nbins),linspace(0,1,nbins));
            xm = 0.5*(xe(1:end-1)+xe(2:end));
            ym = 0.5*(ye(1:end-1)+ye(2:end));
            [X,Y] = meshgrid(xm,ym);
            Z = 1-X-Y;
            
            val(Z<-1e-8) = NaN;
            contourf(v(1,1)*X+v(1,2)*Y+v(1,3)*Z,v(2,1)*X+v(2,2)*Y+v(2,3)*Z,val',"EdgeColor","none")
            axis(ax(ri,ci),"equal")
            xlim(ax(ri,ci),[0 1])
        end
    end
%     for si = 1:total_runs
%         sim_name = ids(si);
% 
%         path_to_sim_folder = sprintf("../data/%s",sim_name);
% 
%         load(sprintf("%s/output_final.mat",path_to_sim_folder),"tracked")
% 
%         a = reshape(tracked.phase_cell_hours,[],4);
%         a = a(:,[1,2,4])';
%         a = a./sum(a,1);
% 
% 
%         %%
%         % scatter3(v(1,:)*a,v(2,:)*a,v(3,:)*a,10,'green','filled')
%         plot(v(1,:)*a,v(2,:)*a,'black',"LineWidth",1)
%     end
end