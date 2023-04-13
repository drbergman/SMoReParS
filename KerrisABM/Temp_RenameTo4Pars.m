% This is a one-off script to add the __CC_12 to all subfolders for
% simulations that ran before we varied dlim. In those runs, it was fixed
% to 12, hence CC_12.

clearvars;

f = dir("data/Binary/AA*");

for i = 1:numel(f)
    s = dir(sprintf("%s/%s/BB*",f(i).folder,f(i).name));
    for j = 1:numel(s)
        I = strfind(s(j).name,'__');
        if length(I)~=1
            error("multiple __ found?")
        end
        new_name = sprintf("%s__CC_12%s",s(j).name(1:I-1),s(j).name(I:end));
        movefile(sprintf("%s/%s",s(j).folder,s(j).name),sprintf("%s/%s",s(j).folder,new_name))
    end
end