
if ~exist("f","var")
    f = dir("./sims/");
end

%%
I = false(size(f));
for i = 1:numel(f)
    I(i) = ~startsWith(f(i).name,".");
end
%%
f = f(I);

%%
I = false(size(f));
for i = 1:numel(f)
    I(i) = startsWith(f(i).name,"2");
end
%%
id = 1;
n = 0;
for i = 1:numel(f)
    if mod(n,10000)==0
        current_sim_folder = sprintf("./sims_%d",id);
        mkdir(current_sim_folder)
        id = id+1;
    end
    source = sprintf("%s/%s",f(i).folder,f(i).name);

    % target = sprintf("/Volumes/MYPASSPORT/UMLaptop/SMoReParS/TumorSpheroidABM/data/sims/%s",f(i).name);
    if exist(source,"dir")
        movefile(source,current_sim_folder);
    end
    n = n+1;
end