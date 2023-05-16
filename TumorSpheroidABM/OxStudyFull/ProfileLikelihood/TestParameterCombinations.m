clearvars;

load("data/Profiles_SMFromData_clean.mat","out")
sm_par_display_names = ["\lambda","\alpha","K","d_{G1/S}","d_{G2/M}","EC50"];

f = figureOnRight;
tiledlayout("flow")

np = size(out,1);
for xi = 1:np
    for yi = (xi+1):np
        for zi = (yi+1):np
            nexttile; hold on;
            scatter3(out{xi}(xi,:),out{xi}(yi,:),out{xi}(zi,:))
            scatter3(out{yi}(xi,:),out{yi}(yi,:),out{yi}(zi,:))
            scatter3(out{zi}(xi,:),out{zi}(yi,:),out{zi}(zi,:))
            xlabel(sm_par_display_names(xi))
            ylabel(sm_par_display_names(yi))
            zlabel(sm_par_display_names(zi))
            view(3)
        end
    end
end
