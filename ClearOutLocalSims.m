if ~exist("f","var")
    f = dir("TumorSpheroidABM/data/sims/");
end

ex_hd_path = "/Volumes/MYPASSPORT/UMLaptop/SMoReParS/TumorSpheroidABM/data/";

if ~exist("ex_hd_dirs","var")
    ex_hd_dirs = dir(ex_hd_path);
end


%% 

for i = 1:numel(f)
    sim_name = f(i).name;

    if startsWith(sim_name,".") || ~exist(sprintf("%s/%s",f(i).folder,sim_name),"dir")
        continue
    end

    sim_files = dir(sprintf("%s/%s",f(i).folder,sim_name));
    for k = numel(sim_files):-1:1
        if startsWith(sim_files(k).name,".")
            sim_files(k) = [];
        end
    end
    
    % check if already on external hard drive
    for j = 1:numel(ex_hd_dirs)
        data_sim_folder = sprintf("%s/%s",ex_hd_dirs(j).folder,ex_hd_dirs(j).name);
        if exist(sprintf("%s/%s",data_sim_folder,sim_name),"dir")
            for k = 1:numel(sim_files)
                if ~exist(sprintf("%s/%s/%s",data_sim_folder,sim_name,sim_files(k).name),"file")
                    movefile(sprintf("%s/%s",sim_files(k).folder,sim_files(k).name),sprintf("%s/%s",data_sim_folder,sim_name))
                end
            end
            rmdir(sprintf("%s/%s",f(i).folder,sim_name),"s")
            break;
        end
    end

end