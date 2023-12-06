clearvars

f1 = dir("test/sims_*");
f2 = dir("TumorSpheroidABM/data/sims_*");
f = [f1;f2];
clear f1 f2

for i = 1:numel(f)
    g = dir(sprintf("%s/%s",f(i).folder,f(i).name));
    for j = 1:numel(g)
        if startsWith(g(j).name,".")
            continue
        end
        destination = sprintf("TumorSpheroidABM/data/sims/%s",g(j).name);
        if exist(g(j).name,"dir")
            continue
        end
        source = sprintf("%s/%s",g(j).folder,g(j).name);
        movefile(source,destination)
    end
end