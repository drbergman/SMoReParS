clearvars;
if ~exist("Data/Binary","dir")
    mkdir("Data/Binary")
end
for i = 1:162
    main(i)
end