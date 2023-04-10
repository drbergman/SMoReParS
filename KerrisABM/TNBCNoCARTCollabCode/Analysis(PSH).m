L = dir;
A = {L.name};
dumbimage = {};
changesinCells = [];
changesinStems = [];
changesinCCR5 = [];
changesinHypo = [];

%Finds all txt files
for i = 1:length(A)
temp_name = A{i};
    if length(temp_name)>= 14
        nor = strcmp(temp_name(1), 'r'); %No rundata
        if strcmp (temp_name(end-3:end), '.txt') && (nor == 0); 
            dumbimage = cat(1,dumbimage,temp_name);
        end
    end
end

%Go through the list and separate into stem and cell files
for j = 1:length(dumbimage)
   %store filename
   filename1 = dumbimage{j};
   %check if it is Stems
   if strcmp(filename1(9:13), 'Stems') 
        %load file
        liststem = dlmread(filename1);
        changesinStems = cat(2,changesinStems, liststem);
        
        
   elseif strcmp (filename1(9:13), 'Cells') 
               %load file
        listcell = dlmread(filename1);
        changesinCells = cat(2,changesinCells, listcell);
       
   elseif strcmp (filename1(9:12), 'CCR5') 
               %load file
        listccr5 = dlmread(filename1);
        changesinCCR5 = cat(2,changesinCCR5, listccr5);
        
   elseif strcmp (filename1(9:12), 'Hypo') 
               %load file
        listhypo = dlmread(filename1);
        changesinHypo = cat(2,changesinHypo, listhypo);
        
   end
    
end

    stem1dataC = changesinCells;
    stem1dataS = changesinStems;
    stem1dataR = changesinCCR5;
    
    changesinHypo(changesinHypo == 0) = NaN;
    stem1dataC(stem1dataC == 0) = NaN;%Don't want zeros to be included in mean
    stem1dataS(stem1dataS == 0) = NaN;
    stem1dataR(stem1dataR == 0) = NaN;

     meanC1 = mean(stem1dataC,2); %Don't take in
    %to account the zeros
errorC1 = std(stem1dataC,0,2);
    meanS1 = mean(stem1dataS,2);
    errorS1 = std(stem1dataS,0,2);
    meanR1 = mean(stem1dataR,2);
    meanlistcell = meanC1;
    meanliststem = meanS1;
    meanlistccr5 = meanR1;
    errorlistcell = errorC1;
    errorliststem = errorS1;

    
figure
plot(meanlistcell')
ylabel('Number of Cells', 'fontsize', 22)
xlabel('Number of Days','fontsize', 22)

figure
plot(changesinHypo')
ylabel('Number of Hypoxic Cells', 'fontsize', 22)
xlabel('Number of Days','fontsize', 22)

figure
plot(meanliststem')
ylabel('Number of Stem Cells','fontsize', 22)
xlabel('Number of Days','fontsize', 22)