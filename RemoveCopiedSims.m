clearvars;

f = dir("/Volumes/MYPASSPORT/UMLaptop/SMoReParS/TumorSpheroidABM/data/");

%% create final folder
for i = numel(f):-1:1
    f_names(i) = string(f(i).name);
end
N = zeros(100,1);
for i = 1:100
    if any(f_names==sprintf("sim_bin_%d",i))
        continue
    end
    dir_name = sprintf("%s/sim_bin_%d",f(1).folder,i);
    mkdir(dir_name)
    temp = dir(dir_name);
    for j = 1:numel(temp)
        N(i) = N(i) + (~startsWith(temp(j).name,"."));
    end
end

%%
for i  = numel(f):-1:1
    if startsWith(f(i).name,".") || startsWith(f(i).name,"sim_bin_")
        continue
    end
    g = dir(sprintf("%s/%s",f(i).folder,f(i).name)); % all the simulation folders
    for j = 1:numel(g)
        if startsWith(g(j).name,".")
            continue
        end
        %% check if local copy exists
        local_copy = sprintf("TumorSpheroidABM/data/sims/%s",g(j).name);
        if exist(local_copy,"dir")
            h = dir(local_copy);
            for k = 1:numel(h)
                destination = sprintf("%s/%s/%s",g(j).folder,g(j).name,h(k).name);
                if startsWith(h(k).name,".") || exist(destination,"file")
                    continue
                end
                source = sprintf("%s/%s",h(k).folder,h(k).name);
                movefile(source,destination)
            end
            rmdir(local_copy,"s")
        end

        %% check if already in a sim_bin
        already_binned = false;
        for k = 1:100
            if exist(sprintf("%s/sim_bin_%d/%s",f(1).folder,k,g(j).name),"dir")
                already_binned = true;
                break
            end
        end

        %% if needed move into sim_bin
        source = sprintf("%s/%s",g(j).folder,g(j).name);
        if already_binned
            rmdir(source)
            continue
        end
        I = find(N < 10000);
        movefile(source,sprintf("%s/sim_bin_%d",f(1).folder,I))
    end
end