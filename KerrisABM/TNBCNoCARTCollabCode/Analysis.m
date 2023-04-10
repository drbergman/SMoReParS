cd Data/AA_10__BB_3
L = dir;
A = {L.name};
dumbimage = {};
changesinCells = [];
changesinStems = [];
changesinCCR5 = [];
changesinHypo = [];
changesinDeath = [];
changesinMacro = [];
changesinKills = [];
changesinKT = [];

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
        
    elseif strcmp (filename1(9:13), 'Macro')
        %load file
        listmacro = dlmread(filename1);
        changesinMacro = cat(2,changesinMacro, listmacro);
        
    elseif strcmp (filename1(9:10), 'KT')
        %load file
        listkt = dlmread(filename1);
        changesinKT = cat(2,changesinKT, listkt);
        
    elseif strcmp (filename1(9:13), 'Death')
        %load file
        listdeath = dlmread(filename1);
        changesinDeath = cat(2,changesinDeath, listdeath);
        
    elseif strcmp (filename1(9:12), 'CCR5')
        %load file
        listccr5 = dlmread(filename1);
        changesinCCR5 = cat(2,changesinCCR5, listccr5);
        
    elseif strcmp (filename1(9:12), 'Hypo')
        %load file
        listhypo = dlmread(filename1);
        changesinHypo = cat(2,changesinHypo, listhypo);
        
    elseif strcmp (filename1(9:13), 'Kills')
        %load file
        listkills = dlmread(filename1);
        changesinKills = cat(2,changesinKills, listkills);
        
    end
    
end

changesinStems(changesinStems == 0) = NaN;
changesinCells(changesinCells == 0) = NaN;
changesinDeath(changesinDeath == 0) = NaN;
changesinHypo(changesinHypo == 0) = NaN;
changesinMacro(changesinMacro == 0) = NaN;
changesinKills(changesinKills == 0) = NaN;
changesinKills(changesinKT == 0) = NaN;



figure
plot(changesinCells)
ylabel('Number of Cells', 'fontsize', 22)
xlabel('Number of Iterations','fontsize', 22)

figure
plot(changesinDeath)
ylabel('Number of Deaths', 'fontsize', 22)
xlabel('Number of Iterations','fontsize', 22)

figure
plot(changesinKills)
ylabel('Number of Kills', 'fontsize', 22)
xlabel('Number of Iterations','fontsize', 22)

figure
plot(changesinKT)
ylabel('Number of CTLs', 'fontsize', 22)
xlabel('Number of Iterations','fontsize', 22)

figure
plot(changesinMacro)
ylabel('Number of Macrophages', 'fontsize', 22)
xlabel('Number of Iterations','fontsize', 22)


figure
plot(changesinHypo)
ylabel('Number of Hypoxic Cells', 'fontsize', 22)
xlabel('Number of Iterations','fontsize', 22)

figure
plot(changesinStems)
ylabel('Number of Stem Cells','fontsize', 22)
xlabel('Number of Iterations','fontsize', 22)

cd ../..