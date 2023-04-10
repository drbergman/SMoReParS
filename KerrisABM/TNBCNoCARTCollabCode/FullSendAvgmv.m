%FullSendAvg:
%A version of FullSend meant to create graphs of multiple runs. ONLY ONE
%TYPE OF GRAPH CAN BE UNCOMMENTED AT A TIME


PERCENT_DEAD = false; %changes num dead to num dead/total cells

BB = [1:6]; %number of runs for the same vars
FOLDERS = 3; %number of runs w different vars (ie folders in BB)

lineSty = {'o','+','*','s','.'};
colorselect = {'r','g','b','c','m'};



cd Data

numvariable = 3;
numtrials = 6;

cellm = cell(numvariable,numtrials);
ktm = cell(numvariable,numtrials);
killsm = cell(numvariable,numtrials);
stemsm = cell(numvariable,numtrials);
ccr5m = cell(numvariable,numtrials);


for a = 1:numvariable
    for b = 1:numtrials
        if a == 1
            aa = 0;
        elseif a == 2
            aa = 10;
        elseif a == numvariable
            aa = 25;
        end
        dirname = strcat('AA_',num2str(aa),'__BB_',num2str(b));
        cd(dirname)
        L = dir;
        A = {L.name};
        for i = 1:length(A)
            temp_name = A{i};
            
            if length(temp_name)<14
                continue
            end
            
            if strcmp (temp_name(9:13), 'Cells')
                listcell = dlmread(temp_name);
                cellm{a,b} = listcell;
                
            elseif strcmp (temp_name(9:10), 'KT')
                listkt = dlmread(temp_name);
                ktm{a,b} = listkt;
                
            elseif strcmp (temp_name(9:13), 'Kills')
                listkills = dlmread(temp_name);
                killsm{a,b} = listkills;
                
            elseif strcmp (temp_name(9:13), 'Stems')
                liststems = dlmread(temp_name);
                stemsm{a,b} = liststems;
                
            elseif strcmp (temp_name(9:12), 'CCR5')
                listccr5 = dlmread(temp_name);
                ccr5m{a,b} = listccr5;
                
                
            end
            
        end
        cd ..
    end
end


%change num dead -> percent dead
if(PERCENT_DEAD)
    for a = AA
        for b = AB/5
            for c = BB
                for ff = 1:length(deathm{a,b,c})
                    if  ~isnan(deathm{a,b,c}(ff))
                        deathm{a,b,c}(ff)=100*(deathm{a,b,c}(ff)/cellm{a,b,c}(ff));
                    end
                end
            end
        end
    end
end

%avging time
newcellm = cell(numvariable,numtrials);
newktm = cell(numvariable,numtrials);
newkillsm = cell(numvariable,numtrials);
newstemsm = cell(numvariable,numtrials);
newccr5m = cell(numvariable,numtrials);

fullcellm = cell(numvariable,numtrials);
fullktm = cell(numvariable,numtrials);
fullkillsm = cell(numvariable,numtrials);
fullstemsm = cell(numvariable,numtrials);
fullccr5m = cell(numvariable,numtrials);

for a = 1:numvariable
        temp = zeros(300,5);
        z = zeros(300,numtrials);
        x = zeros(300,numtrials);
        y = zeros(300,numtrials);
        q = zeros(300,numtrials);
        w = zeros(300,numtrials);
        for b = numtrials
            C1 = cellm{a,b};
            z(:,b) = C1;
            
            C2 = ktm{a,b};
            x(:,b) = C2;
            
            C3 = killsm{a,b};
            y(:,b) = C3;
            
            C4 = stemsm{a,b};
            q(:,b) = C4;
            
            C5 = ccr5m{a,b};
            w(:,b) = C5;
            
            for d = 1:length(C1)
                temp(d,1) = temp(d,1) + C1(d);
            end
            
            for d = 1:length(C2)
                temp(d,2) = temp(d,2) + C2(d);
            end
            
            for d = 1:length(C3)
                temp(d,3) = temp(d,3) + C3(d);
            end
            
            for d = 1:length(C4)
                temp(d,4) = temp(d,4) + C4(d);
            end
            
            for d = 1:length(C5)
                temp(d,5) = temp(d,5) + C5(d);
            end
        end
        newcellm(a,b) = {temp(:,1)./length(BB)};
        fullcellm(a,b) = {z};
        
        newktm(a,b) = {temp(:,2)./length(BB)};
        fullktm(a,b) = {x};
        
        newkillsm(a,b) = {temp(:,3)./length(BB)};
        fullkillsm(a,b) = {y};
        
        newstemsm(a,b) = {temp(:,4)./length(BB)};
        fullstemsm(a,b) = {q};
        
        newccr5m(a,b) = {temp(:,5)./length(BB)};
        fullccr5m(a,b) = {w};
end

newcellm(a,b);
newktm(a,b);
newkillsm(a,b);
newstemsm(a,b);
newccr5m(a,b);

%calc confidence
N = 6;
zVal = 1.96; % 95%
for a = 1:numvariable
    for b = numtrials
        devsC = zeros(300,1);
        confValsC = zeros(300,1);
        
        devsKT = zeros(300,1);
        confValsKT = zeros(300,1);
        
        devsK = zeros(300,1);
        confValsK = zeros(300,1);
        
        devsS = zeros(300,1);
        confValsS = zeros(300,1);
        
        devs5 = zeros(300,1);
        confVals5 = zeros(300,1);
        
        for ii = 1:300
            devsC(ii) = std(fullcellm{a,b}(ii,:));
            confValsC(ii) = zVal*(devsC(ii)/sqrt(N));
            
            devsKT(ii) = std(fullktm{a,b}(ii,:));
            confValsKT(ii) = zVal*(devsKT(ii)/sqrt(N));
            
            devsK(ii) = std(fullkillsm{a,b}(ii,:));
            confValsK(ii) = zVal*(devsK(ii)/sqrt(N));
            
            devsS(ii) = std(fullstemsm{a,b}(ii,:));
            confValsS(ii) = zVal*(devsS(ii)/sqrt(N));
            
            devs5(ii) = std(fullccr5m{a,b}(ii,:));
            confVals5(ii) = zVal*(devs5(ii)/sqrt(N));
        end
    end
end

%%UNCOMMENT FOR PLOTTING NUM CELLS

figure
ylabel('Number of Cells', 'fontsize', 22)
xlabel('Number of Iterations','fontsize', 22)
hold on
for c = 1:numvariable
    for d = numtrials
        if c == 1
            cc = 0;
        elseif c == 2
            cc = 10;
        elseif c == numvariable
            cc = 25;
        end
            plot(cellm{c,d},strcat(lineSty{c},colorselect{c}),'DisplayName',strcat('# of CTLs: ',num2str(cc)));
            %plot(bounds(:,1),strcat('^',colorselect{FF}),'HandleVisibility','off','MarkerSize',3);
            %plot(bounds(:,2),strcat('v',colorselect{FF}),'HandleVisibility','off','MarkerSize',3);
            errorbar(cellm{c,d}, confValsC,colorselect{c}, 'LineStyle','none','HandleVisibility','off');
    end
end
hold off
legend('Location','northwest')
%         saveas(F,strcat("CellsAvg_AA_",num2str(a),"_AB_",num2str(b*5),'.png'));

figure
ylabel('Number of CTLs', 'fontsize', 22)
xlabel('Number of Iterations','fontsize', 22)
hold on
for c = 1:numvariable
    for d = numtrials
        if c == 1
            cc = 0;
        elseif c == 2
            cc = 10;
        elseif c == numvariable
            cc = 25;
        end
            plot(ktm{c,d},strcat(lineSty{c},colorselect{c}),'DisplayName',strcat('# of CTLs: ',num2str(cc)));
            %plot(bounds(:,1),strcat('^',colorselect{FF}),'HandleVisibility','off','MarkerSize',3);
            %plot(bounds(:,2),strcat('v',colorselect{FF}),'HandleVisibility','off','MarkerSize',3);
            errorbar(ktm{c,d}, confValsKT,colorselect{c}, 'LineStyle','none','HandleVisibility','off');
    end
end
hold off
legend('Location','northwest')
%         saveas(F,strcat("CellsAvg_AA_",num2str(a),"_AB_",num2str(b*5),'.png'));

figure
ylabel('Number of CTL Kills', 'fontsize', 22)
xlabel('Number of Iterations','fontsize', 22)
hold on
for c = 1:numvariable
    for d = numtrials
        if c == 1
            cc = 0;
        elseif c == 2
            cc = 10;
        elseif c == numvariable
            cc = 25;
        end
            plot(killsm{c,d},strcat(lineSty{c},colorselect{c}),'DisplayName',strcat('# of CTLs: ',num2str(cc)));
            %plot(bounds(:,1),strcat('^',colorselect{FF}),'HandleVisibility','off','MarkerSize',3);
            %plot(bounds(:,2),strcat('v',colorselect{FF}),'HandleVisibility','off','MarkerSize',3);
            errorbar(killsm{c,d}, confValsK,colorselect{c}, 'LineStyle','none','HandleVisibility','off');
    end
end
hold off
legend('Location','northwest')
%         saveas(F,strcat("CellsAvg_AA_",num2str(a),"_AB_",num2str(b*5),'.png'));

figure
ylabel('Number of Stem Cells', 'fontsize', 22)
xlabel('Number of Iterations','fontsize', 22)
hold on
for c = 1:numvariable
    for d = numtrials
        if c == 1
            cc = 0;
        elseif c == 2
            cc = 10;
        elseif c == numvariable
            cc = 25;
        end
            plot(stemsm{c,d},strcat(lineSty{c},colorselect{c}),'DisplayName',strcat('# of CTLs: ',num2str(cc)));
            %plot(bounds(:,1),strcat('^',colorselect{FF}),'HandleVisibility','off','MarkerSize',3);
            %plot(bounds(:,2),strcat('v',colorselect{FF}),'HandleVisibility','off','MarkerSize',3);
            errorbar(stemsm{c,d}, confValsS,colorselect{c}, 'LineStyle','none','HandleVisibility','off');
    end
end
hold off
legend('Location','northwest')
%         saveas(F,strcat("CellsAvg_AA_",num2str(a),"_AB_",num2str(b*5),'.png'));

figure
hold
ylabel('Number of CCR5+ Cells', 'fontsize', 22)
xlabel('Number of Iterations','fontsize', 22)
hold on
for c = 1:numvariable
    for d = numtrials
        if c == 1
            cc = 0;
        elseif c == 2
            cc = 10;
        elseif c == numvariable
            cc = 25;
        end
            plot(ccr5m{c,d},strcat(lineSty{c},colorselect{c}),'DisplayName',strcat('# of CTLs: ',num2str(cc)));
            %plot(bounds(:,1),strcat('^',colorselect{FF}),'HandleVisibility','off','MarkerSize',3);
            %plot(bounds(:,2),strcat('v',colorselect{FF}),'HandleVisibility','off','MarkerSize',3);
            errorbar(ccr5m{c,d}, confVals5,colorselect{c}, 'LineStyle','none','HandleVisibility','off');
    end
end
hold off
legend('Location','northwest')
%         saveas(F,strcat("CellsAvg_AA_",num2str(a),"_AB_",num2str(b*5),'.png'));


disp('DONE!!!')
cd ..
% save("fullDat.mat","fullDat");