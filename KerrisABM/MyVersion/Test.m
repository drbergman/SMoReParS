clearvars;
if ~exist("Data/Binary","dir")
    mkdir("Data/Binary")
end
for i = 28
    main(i)
end