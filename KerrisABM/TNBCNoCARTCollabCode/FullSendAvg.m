%FullSendAvg:
%A version of FullSend meant to create graphs of multiple runs. ONLY ONE
%TYPE OF GRAPH CAN BE UNCOMMENTED AT A TIME


PERCENT_HYP = true; %chnages hypoxic to hypoxic/total cells
PERCENT_DEAD = false; %changes num dead to num dead/total cells
JET_FIG = false;
SEN_DEATH = false; %changes death to death-(hypoxic/10) giving us deaths that are not due to hypoxia

BB = [1:6]; %number of runs for the same vars
F = figure;
G = figure;
H = figure;
I = figure;
cd Data/BB
FOLDERS = 3; %number of runs w different vars (ie folders in BB)


for FF = 3:FOLDERS+2
    DIR = struct2cell(dir);
    curF = DIR{1,FF};
    cd(curF);
    AA = [str2num(curF(1))];
    if length(curF)==9
        AB = [str2num(curF(end-6:end-5))];
    else
        AB = [str2num(curF(end-5:end-5))];
    end
    varyStem = true;
    
    cellm = cell(5,4,6);
    hypm = cell(5,4,6);
    deathm = cell(5,4,6);
    macrom = cell(5,4);
    stemm = cell(5,4);
    
    
    lineSty = {'o','+','*','s','.'};
    colorselect = {'r','g','b','c','m'};
    
    %gather data
    for c = BB
        for a = AA
            for b = AB
                b= b/5;
                dirname = strcat('AA_',num2str(a),'__AB_',num2str(b*5),'__BB_',num2str(c));
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
                        listcell(listcell == 0) = NaN;
                        cellm{a,b,c} = listcell;
                        
                    elseif strcmp (temp_name(9:13), 'Macro')
                        listmacro = dlmread(temp_name);
                        listmacro(listmacro == 0) = NaN;
                        macrom{a,b} = listmacro;
                        
                    elseif strcmp (temp_name(9:13), 'Death')
                        listdead = dlmread(temp_name);
                        listdead(listdead == 0) = NaN;
                        deathm{a,b,c} = listdead;
                        
                    elseif strcmp (temp_name(9:12), 'Hypo')
                        listhypo = dlmread(temp_name);
                        listhypo(listhypo == 0) = NaN;
                        hypm{a,b,c} = listhypo;
                        
                    elseif strcmp (temp_name(9:12), 'Stem')
                        liststem = dlmread(temp_name);
                        liststem(liststem == 0) = NaN;
                        stemm{a,b} = liststem;
                        
                    end
                    
                end
                cd ..
            end
        end
    end
    cd ../../..
    cd Plots
    
    if SEN_DEATH
        for a = AA
            for b = AB/5
                for c = BB
                    for ff = 1:length(deathm{a,b,c})
                        if  ~isnan(deathm{a,b,c}(ff))
                            deathm{a,b,c}(ff) = deathm{a,b,c}(ff) - (hypm{a,b,c}(ff)/10);
                        end
                    end
                end
            end
        end
    end
    
    %change num hyp -> percent hyp
    if(PERCENT_HYP)
        for a = AA
            for b = AB/5
                for c = BB
                    for ff = 1:length(hypm{a,b,c})
                        if  ~isnan(hypm{a,b,c}(ff))
                            hypm{a,b,c}(ff)=100*(hypm{a,b,c}(ff)/cellm{a,b,c}(ff));
                        end
                    end
                end
            end
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
    newcellm = cell(5,4);
    fullDat = cell(5,4);
    
    newhypm = cell(5,4);
    fullhyp = cell(5,4);
    
    newdeathm = cell(5,4);
    fulldeath = cell(5,4);
    
    for a = AA
        for b = AB/5
            temp = zeros(300,3);
            z = zeros(300,6);
            x = zeros(300,6);
            y = zeros(300,6);
            for c = BB
                Cell = cellm{a,b,c};
                z(:,c) = Cell;
                
                C2 = hypm{a,b,c};
                x(:,c) = C2;
                
                C3 = deathm{a,b,c};
                y(:,c) = C3;
                
                for d = 1:length(Cell)
                    temp(d,1) = temp(d,1) + Cell(d);
                end
                
                for d = 1:length(C2)
                    temp(d,2) = temp(d,2) + C2(d);
                end
                
                for d = 1:length(C3)
                    temp(d,3) = temp(d,3) + C3(d);
                end
            end
            newcellm(a,b) = {temp(:,1)./length(BB)};
            fullDat(a,b) = {z};
            
            newhypm(a,b) = {temp(:,2)./length(BB)};
            fullhyp(a,b) = {x};
            
            newdeathm(a,b) = {temp(:,3)./length(BB)};
            fulldeath(a,b) = {y};
        end
    end
    cellm(a,b) = newcellm(a,b);
    hypm(a,b) = newhypm(a,b);
    deathm(a,b) = newdeathm(a,b);
    
    %calc confidence
    N = 6;
    zVal = 1.96; % 95%
    for a = AA
        for b = AB./5
            devs = zeros(300,1);
            confVals = zeros(300,1);
            
            devsH = zeros(300,1);
            confValsH = zeros(300,1);
            
            devsD = zeros(300,1);
            confValsD = zeros(300,1);
            for ii = 1:300
                devs(ii) = std(fullDat{a,b}(ii,:));
                confVals(ii) = zVal*(devs(ii)/sqrt(N));
                
                devsH(ii) = std(fullhyp{a,b}(ii,:));
                confValsH(ii) = zVal*(devsH(ii)/sqrt(N));
                
                devsD(ii) = std(fulldeath{a,b}(ii,:));
                confValsD(ii) = zVal*(devsD(ii)/sqrt(N));
            end
        end
    end
    
    
    if JET_FIG
        jetdat = cell(5,4);
        jetdaty = cell(5,4);
        for a = AA
            for b = AB
                row = [];
                rowy = [];
                curr = 20;
                for jj = 1:length(stemm{a,b})
                    if stemm{a,b}(jj)>curr
                        row(end+1) = jj;
                        rowy(end+1) = cellm{a,b}(jj);
                        curr = stemm{a,b}(jj);
                    end
                end
                jetdat{a,b}=row;
                jetdaty{a,b} = rowy;
            end
        end
    end
    
    %%UNCOMMENT FOR PLOTTING NUM CELLS
    
        ylabel('Number of Cells', 'fontsize', 22)
        xlabel('Number of Iterations','fontsize', 22)
        hold on
        for c = AA
            for d = AB/5
                if JET_FIG
                    if varyStem
                        plot(cellm{c,d},strcat(lineSty{c},colorselect{c}),'DisplayName',strcat('Stem ',num2str(c/100),' & VEGF ',num2str(d*5)),'MarkerSize',4);
                        plot(jetdat{c,d},jetdaty{c,d},strcat("d",colorselect{c}),'MarkerSize',8,'HandleVisibility','off','MarkerFaceColor',colorselect{c});
                    else
                        plot(cellm{c,d},strcat(lineSty{d},colorselect{d}),'DisplayName',strcat('Stem ',num2str(c/100),' & VEGF ',num2str(d*5)),'MarkerSize',4);
                        plot(jetdat{c,d},jetdaty{c,d},strcat("d",colorselect{d}),'MarkerSize',8,'HandleVisibility','off','MarkerFaceColor',colorselect{d});
                    end
                else
                    plot(cellm{c,d},strcat(lineSty{d},colorselect{FF}),'DisplayName',strcat('Stem ',num2str(c/100),' & VEGF ',num2str(d*5)));
                    %plot(bounds(:,1),strcat('^',colorselect{FF}),'HandleVisibility','off','MarkerSize',3);
                    %plot(bounds(:,2),strcat('v',colorselect{FF}),'HandleVisibility','off','MarkerSize',3);
                    errorbar(cellm{c,d}, confVals,colorselect{FF}, 'LineStyle','none','HandleVisibility','off');
                end
            end
        end
        legend('Location','northwest')
%         saveas(F,strcat("CellsAvg_AA_",num2str(a),"_AB_",num2str(b*5),'.png'));
    

%%UNCOMMENT FOR PLOTTING HYPOXIC

%         if PERCENT_HYP
%             ylabel('% of Hypoxic Cells', 'fontsize', 22)
%         else
%             ylabel('Number of Hypoxic Cells', 'fontsize', 22)
%         end
%         xlabel('Number of Iterations','fontsize', 22)
%         hold on
%         for c = AA
%             for d = AB/5
%                 plot(hypm{c,d},strcat(lineSty{d},colorselect{FF}),'DisplayName',strcat('Stem ',num2str(c),' & VEGF ',num2str(d*5)))
%                 errorbar(hypm{c,d}, confValsH,colorselect{FF}, 'LineStyle','none','HandleVisibility','off');
%             end
%         end
%         legend('Location','northwest')
    

%%UNCOMMENT FOR PLOTTING DEAD

%     if PERCENT_DEAD
%         ylabel('% of Cell Deaths', 'fontsize', 22)
%     else
%         ylabel('Number of Cell Deaths', 'fontsize', 22)
%     end
%     xlabel('Number of Iterations','fontsize', 22)
%     hold on
%     for c = AA
%         for d = AB/5
%             plot(deathm{c,d},strcat(lineSty{d},colorselect{FF}),'DisplayName',strcat('Stem ',num2str(c),' & VEGF ',num2str(d*5)))
%             errorbar(newdeathm{c,d}, confValsD,colorselect{FF}, 'LineStyle','none','HandleVisibility','off');
%         end
%     end
%     legend('Location','northwest')
    

%%UNCOMMENT FOR PLOTTING NUM MACROPHAGES

    % ylabel('Number of Macrophages', 'fontsize', 22)
    % xlabel('Number of Iterations','fontsize', 22)
    % hold on
    % for c = 1:5
    %     for d = 1:4
    %         plot(macrom{c,d},lineSty{d},'DisplayName',strcat('Stem ',num2str(c),' & VEGF ',num2str(d*5)))
    %     end
    % end
    % legend('Location','northwest')
    %
    % ylabel('Number of Stem Cells', 'fontsize', 22)
    % xlabel('Number of Iterations','fontsize', 22)
    % hold on
    % for c = 1:5
    %     for d = 1:4
    %         plot(stemm{c,d},lineSty{d},'DisplayName',strcat('Stem ',num2str(c),' & VEGF ',num2str(d*5)))
    %     end
    % end
    % legend('Location','northwest')
    
    disp('DONE!!!')
    cd ../Data/BB
end
save("fullDat.mat","fullDat");