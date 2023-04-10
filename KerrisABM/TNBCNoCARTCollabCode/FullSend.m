%cd Data
cd StemOnlyData

numvariable = 1;
numtrials = 2;

cellm = cell(numvariable,numtrials);
ktm = cell(numvariable,numtrials);
killsm = cell(numvariable,numtrials);
stemsm = cell(numvariable,numtrials);


for a = numvariable
    for b = 2:numtrials
        if a == 1
            aa = 0;
        elseif a == 2
            aa = 10;
        elseif a == numvariable
            aa = 25;
        end
        aa=1%rm after
        dirname = strcat('AA_',num2str(aa),'__BB_',num2str(b));
        cd(dirname)
        L = dir;
        A = {L.name};
        for i = 1:length(A)
            temp_name = A{i};
            
            if length(temp_name)<14
                continue
            end
            
            if strcmp (temp_name(9:13), 'Cells')%was 12 -> 13
                listcell = dlmread(temp_name);
                listcell(listcell == 0) = NaN;
                cellm{a,b} = listcell;
                
            elseif strcmp (temp_name(9:10), 'KT')
                listkt = dlmread(temp_name);
                listkt(listkt == 0) = NaN;
                ktm{a,b} = listkt;
                
            elseif strcmp (temp_name(9:13), 'Kills')
                listkills = dlmread(temp_name);
                listkills(listkills == 0) = NaN;
                killsm{a,b} = listkills;
                
            elseif strcmp (temp_name(9:13), 'Stems')
                liststems = dlmread(temp_name);
                liststems(liststems == 0) = NaN;
                stemsm{a,b} = liststems;
                
            end
            
        end
        cd ..
    end
end

figure
ylabel('Number of Cells', 'fontsize', 22)
xlabel('Number of Iterations','fontsize', 22)
hold on
for c = numvariable
    for d = 1:numtrials
        if c == 1
            bb = 0;
        elseif a == 2
            bb = 10;
        elseif a == numvariable
            bb = 25;
        end
        plot(cellm{c,d},'DisplayName',strcat('# of CTLs: ',num2str(bb),' Trial #',num2str(d)))
    end
end
hold off
legend('Location','northwest')

figure
ylabel('Number of CTLs', 'fontsize', 22)
xlabel('Number of Iterations','fontsize', 22)
hold on
for c = numvariable
    for d = 1:numtrials
        if c == 1
            bb = 0;
        elseif a == 2
            bb = 10;
        elseif a == numvariable
            bb = 25;
        end
        plot(ktm{c,d},'DisplayName',strcat('# of CTLs: ',num2str(bb),' Trial #',num2str(d)))
    end
end
hold off
legend('Location','northwest')

figure
ylabel('Number of CTL Kills', 'fontsize', 22)
xlabel('Number of Iterations','fontsize', 22)
hold on
for c = numvariable
    for d = 1:numtrials
        if c == 1
            bb = 0;
        elseif a == 2
            bb = 10;
        elseif a == numvariable
            bb = 25;
        end
        plot(killsm{c,d},'DisplayName',strcat('# of CTLs: ',num2str(bb),' Trial #',num2str(d)))
    end
end
hold off
legend('Location','northwest')

figure
ylabel('Number of Stem Cells', 'fontsize', 22)
xlabel('Number of Iterations','fontsize', 22)
hold on
for c = numvariable
    for d = 1:numtrials
        if c == 1
            bb = 0;
        elseif a == 2
            bb = 10;
        elseif a == numvariable
            bb = 25;
        end
        plot(stemsm{c,d},'DisplayName',strcat('# of CTLs: ',num2str(bb),' Trial #',num2str(d)))
    end
end
hold off
legend('Location','northwest')

cd ..