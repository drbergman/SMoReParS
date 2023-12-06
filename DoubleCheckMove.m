% clearvars;

if ~exist("f","var")
f = dir("TumorSpheroidABM/data/sims/");
end

%%
for i = 131285:numel(f)
    if mod(i,1000)==0
        fprintf("Finished %d.\n",i)
    end
    if startsWith(f(i).name,".")
        continue
    end
    destination = sprintf("/Volumes/MYPASSPORT/UMLaptop/SMoReParS/TumorSpheroidABM/data/all_sims/%s",f(i).name);
    if exist(destination,"dir")
        continue
    end
    fprintf("  Moving %s.\n",f(i).name)
    source = sprintf("%s/%s",f(i).folder,f(i).name);
    try
        movefile(source,destination)
    catch % sometimes it doesn't work?
        mkdir(destination)
        h = dir(source);
        for j = 1:numel(h)
            if startsWith(h(j).name,".")
                continue
            end
            movefile(sprintf("%s/%s",h(j).folder,h(j).name),destination)
        end
    end
end

