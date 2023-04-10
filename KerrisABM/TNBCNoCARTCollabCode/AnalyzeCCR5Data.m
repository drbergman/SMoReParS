%AnalyzeCCR5Data
%This one uses regular expressions

clear all;
close all;

%Read in files

L = dir;
A = {L.name};
dumbimage = {};
changesinCells = [];
changesinStems = [];
changesinCCR5 = [];
stemrun = [];
stemparam = [];
cellrun = [];
cellparam = [];
CCR5run = [];
CCR5param = [];
bboxmax = {};

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
        stemrun(size(changesinStems,2)) = str2double(filename1(15));
        %regular expression
         express = 'd.*tx';
        Match =regexp(filename1,express,'match');
        Mstr = Match{1};
        stemparam(length(stemrun)) = str2double(Mstr(2:end-3)); %Need if is not a fraction
        
   elseif strcmp (filename1(9:13), 'Cells') 
               %load file
        listcell = dlmread(filename1);
        changesinCells = cat(2,changesinCells, listcell);
        cellrun(size(changesinCells,2)) = str2double(filename1(15));
        %regular expression
        express = 'd.*tx';
        Match =regexp(filename1,express,'match');
        Mstr = Match{1};
        cellparam(length(cellrun)) = str2double(Mstr(2:end-3)); %Need if is not a fraction
        
   elseif strcmp (filename1(9:12), 'CCR5') 
               %load file
        listccr5 = dlmread(filename1);
        changesinCCR5 = cat(2,changesinCCR5, listccr5);
        CCR5run(size(changesinCCR5,2)) = str2double(filename1(15));
        %regular expression
        express = 'd.*tx';
        Match =regexp(filename1,express,'match');
        Mstr = Match{1};
        CCR5param(length(CCR5run)) = str2double(Mstr(2:end-3)); %Need if is not a fraction
   end
    
end

%Remove bad data
% stemparam(1:2) = [];
% cellparam(1:2) = [];
% stemrun(1:2) = [];
% cellrun(1:2) = [];
% 
% changesinCells(:,1:2) = [];
% changesinStems(:,1:2) = [];

%unique params
theparams = unique(cellparam);

  %Find first run data
    firstrun = cellrun == 1;
    firstrundataC = changesinCells(:,firstrun);
    firstrundataS = changesinStems(:,firstrun);
    firstrundataCR = changesinCCR5(:,firstrun);


%Data has been collected - Now graph it!
figure
semilogy(changesinCells,'LineWidth',2)
titlehand =title('All Cell Data')
set(titlehand,'FontSize',24, 'Fontname', 'Arial')
%     title('Mean Cell Data with Changing Stem Cell Symmetric Division Rate')
    lparams = num2str((theparams./100)');
    [legh,objh,outh,outm] = legend(lparams, 'Location', 'SouthEast');
    %title('Symmetric Division Rate'
    axis([0 1000 0 max(max(firstrundataC))])
    v = get(legh,'title');
    set(v,'string','Seeding Rate','FontSize',14, 'Fontname', 'Arial');
    yhand = ylabel('Number of Total Cells');
    xhand = xlabel('Number of Days');
    set(yhand,'FontSize',22, 'Fontname', 'Arial')
    set(xhand,'FontSize',22, 'Fontname', 'Arial')
    set(gca,'FontSize',16, 'Fontname', 'Arial')

figure
semilogy(changesinStems,'LineWidth',2)
titlehand =title('All Stem Data')
set(titlehand,'FontSize',24, 'Fontname', 'Arial')
%     title('Mean Cell Data with Changing Stem Cell Symmetric Division Rate')
    lparams = num2str((theparams./100)');
    [legh,objh,outh,outm] = legend(lparams, 'Location', 'SouthEast');
    %title('Symmetric Division Rate'
    axis([0 1000 0 max(max(firstrundataS))])
    v = get(legh,'title');
    set(v,'string','Seeding Rate','FontSize',14, 'Fontname', 'Arial');
    yhand = ylabel('Number of Total Stems');
    xhand = xlabel('Number of Days');
    set(yhand,'FontSize',22, 'Fontname', 'Arial')
    set(xhand,'FontSize',22, 'Fontname', 'Arial')
    set(gca,'FontSize',16, 'Fontname', 'Arial')

figure
[axes_handles] = semilogy(changesinCCR5, 'LineWidth',2);
titlehand = title('All CCR5 Data');
set(titlehand,'FontSize',24, 'Fontname', 'Arial')
%     title('Mean Cell Data with Changing Stem Cell Symmetric Division Rate')
    lparams = ['Control'; 'Hypoxia'; 'Maraviroc'];%num2str((theparams./100)');
    [legh,objh,outh,outm] = legend(lparams, 'Location', 'SouthEast');
    %title('Symmetric Division Rate'
    axis([0 1000 0 ])
    v = get(legh,'title');
    set(v,'string','Seeding Rate','FontSize',14, 'Fontname', 'Arial');
    yhand = ylabel('Number of Total CCR5');
    xhand = xlabel('Number of Days');
    set(yhand,'FontSize',22, 'Fontname', 'Arial')
    set(xhand,'FontSize',22, 'Fontname', 'Arial')
    set(gca,'FontSize',16, 'Fontname', 'Arial')
    
  




    %Find first run data
    firstrun = cellrun == 2;
    firstrundataC = changesinCells(:,firstrun);
    firstrundataS = changesinStems(:,firstrun);

    %Plot first data
    figure
    plot(firstrundataC)
    figtitle = ['First Run Cell Data with Changing Seeding Rate and No Migration'];
    title(figtitle)
%     title('Mean Cell Data with Changing Stem Cell Symmetric Division Rate')
    lparams = num2str((theparams./100)');
    legend(lparams, 'Location', 'NorthWest')
    %title('Symmetric Division Rate')
    axis([0 1000 0 max(max(firstrundataC))])
    ylabel('Number of Total Cells')
    xlabel('Number of Days')

    figure
    plot(firstrundataS)
    title('First Run Stem Data with Changing Seeding Rate and No Migration')
    lparams = num2str((theparams./100)');
    legend(lparams, 'Location', 'NorthWest')
    %title('Symmetric Division Rate')
    axis([0 1000 0 max(max(firstrundataS))])
    ylabel('Number of Stem Cells')
    xlabel('Number of Days')
    
    
    %preallocate
    meanlistcell = zeros(length(theparams), length(listcell)); 
    meanliststem = zeros(length(theparams), length(liststem)); 
    errorlistcell = zeros(length(theparams), length(listcell)); 
    errorliststem = zeros(length(theparams), length(liststem)); 
    
for cycle = 1:length(theparams)
    param = theparams(cycle);
    %Plot stem = .01 data
    stem1 = cellparam == param;
    stem1dataC = changesinCells(:,stem1);
    stem1dataS = changesinStems(:,stem1);
    stem1dataC(stem1dataC == 0) = NaN;%Don't want zeros to be included in mean
    stem1dataS(stem1dataS == 0) = NaN;
    meanC1 = mean(stem1dataC,2); %Don't take into account the zeros
    errorC1 = std(stem1dataC,0,2);
    meanS1 = mean(stem1dataS,2);
    errorS1 = std(stem1dataS,0,2);
    meanlistcell(cycle,:) = meanC1;
    meanliststem(cycle,:) = meanS1;
    errorlistcell(cycle,:) = errorC1;
    errorliststem(cycle,:) = errorS1;
end
dlmwrite('231HMeanCellList.txt',meanlistcell);
dlmwrite('231HMeanStemList.txt',meanliststem);
dlmwrite('231HMeanCCR5List.txt',meanlistccr5);     

%Plot stem = .03 data
figure
plot(meanlistcell')
title('Division Limit of 6 - No Migration','fontsize', 24)
lparams = num2str((theparams./100)');
[legh,objh,outh,outm] = legend(lparams, 'Location', 'NorthWest');
set(objh,'linewidth',2);
v = get(legh,'title');
set(v,'string','Seeding Rate');
axis([0 1000 0 max(max(meanlistcell))])
ylabel('Mean Number of Cells', 'fontsize', 22)
xlabel('Number of Days','fontsize', 22)

figure
plot(meanliststem')
title('Stem Data of Stem Cell Seeding Rate - Migration','fontsize', 24)
lparams = num2str((theparams./100)');
[legh,objh,outh,outm] = legend(lparams, 'Location', 'NorthWest');
v = get(legh,'title');
set(v,'string','Seeding Rate');
set(objh,'linewidth',2);
axis([0 1000 0 max(max(meanliststem))])
ylabel('Mean Number of Stem Cells','fontsize', 22)
xlabel('Number of Days','fontsize', 22)


% %Plot stem = .01 data
% stem5 = cellparam == 5;
% stem5dataC = changesinCells(:,stem5);
% stem5dataS = changesinStems(:,stem5);
% 
% %Plot stem = .01 data
% stem7 = cellparam ==7;
% stem7dataC = changesinCells(:,stem7);
% stem7dataS = changesinStems(:,stem7);
% 
% %Plot first data
% figure
% plot(stem1dataC)
% title('Stem1 Cell Data')
% 
% 
% 
% %Plot stem = .09 data
% stem9 = cellparam ==9;
% stem9dataC = changesinCells(:,stem9);
% stem9dataS = changesinStems(:,stem9);
% 
% %Plot first data
% figure
% plot(stem9dataC)
% title('Stem9 Cell Data')
% 
% figure
% plot(stem9dataS)
% title('Stem9 Stem Data')
% 
% %Plot means and error bars
% errorbar(y,e,'xr')
