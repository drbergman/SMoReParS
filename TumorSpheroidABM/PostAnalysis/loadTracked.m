function tracked = loadTracked(ids)

for i = numel(ids):-1:1
     temp = load(sprintf("../data/%s/output_final.mat",ids(i)),"tracked");
     tracked(i) = temp.tracked;
end

tracked = reshape(tracked,size(ids));