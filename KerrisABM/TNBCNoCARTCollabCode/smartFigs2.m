%Figure plotting script for experiments run during BSRI 2021.

%Plotting params and options===============================================

%param dim 1
aa = [25,50,75,98];
%param dim 2
ab = [150];
%numer of trials (runs)
bb = 6;

param1 = 'Antigen % '; %antigen expression
param2 = 'CAR-T delay ';

%Options:
%set options to 1 if disired, 0 otherwise.
itters = 300; %number of tics in the simmulation
exType = 'Binary'; %Experiment type: Gradated or Binary
showErrorBars = 1; %only applies to multiple runs
plotRuns = 0; %plots every run for a point in the param space in the same fig
plotAvg = 0; %averages across runs, plots every point in param space in the same fig
saveImgs = 0;
runStats = 1;
ANOVAsampleTime = 299;
%color scheme for figs
colors = [238,204,102;102,153,204;153,68,85;0,68,136;0,0,0]./255;

%colors = hsv(4);

%add data files that should not be plotted to this array (a unique substring of the file name will do):
ignore = ["State","Macro","Vasc","HypCell","Antigen"];

labels = {'Number of CAR-T cells','Number of cancer cells with CCR5 expression','Number of cancer cells','Total number of cancer cell deaths','Number of hypoxic cancer cells','Number of cancer cells killed by CAR-T','Number of cancer proliferation','Number of cancer stem cells','Vascular density'};


%==========================================================================
%Data is stored in a 5-D array
% dim 1 = aa values (len = 4)
% dim 2 = ab values (len = 1)
% dim 3 = bb values (runs) (len = 7)
% dim 4 = all figure types to plot (len = 8)
% dim 5 = itterations of the runs (len = 300)
% lengths of each dim subject to change

%Loading data
%cd to proper directory
cd Data
cd(exType);
exFiles = strcat('AA_',num2str(aa(1)),'__AB_',num2str(ab(1)));
cd(exFiles);
runFiles = strcat('BB_',num2str(1));
cd(runFiles);
cd('fig_data');
datafiles = dir(fullfile('*.txt'));  

%finding files to plot and verifying they exist 
files2plot = [];
for f = 1:length(datafiles)
    pass = 0;
    for ig = ignore
        if contains(datafiles(f).name, ig);
            pass = 1;
        end
    end
    if pass == 1
        continue
    end
    try
        dat = dlmread(datafiles(f).name); 
        files2plot = [files2plot, convertCharsToStrings(datafiles(k).name)];
    catch
        disp(strcat('datafile_:_',datafiles(f).name,' is unreadable'));
    end
end
%definind and preallocating data matrix
data = zeros(length(aa),length(ab),bb,length(files2plot),itters);
cd ../../..

%looping over param space
for a_a = 1:length(aa)
    for a_b = 1:length(ab)
        %looping over each trial
        
        exFiles = strcat('AA_',num2str(aa(a_a)),'__AB_',num2str(ab(a_b)));
        cd(exFiles);
        %looping over each run
        for run = 1:bb
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %cd to proper subdirectory and getting desired data files
            runFiles = strcat('BB_',num2str(run));
            cd(runFiles);
            cd('fig_data');
            datafiles = dir(fullfile('*.txt'));    
            files2plot = [];
            for f = 1:length(datafiles)
                %checking for ignored files
                pass = 0;
                for ig = ignore
                    if contains(datafiles(f).name, ig);
                        pass = 1;
                    end
                end
                if pass == 1
                    continue
                end
                %getting data file names
                try
                    dat = dlmread(datafiles(f).name); 
                    files2plot = [files2plot, convertCharsToStrings(datafiles(f).name)];
                catch
                    disp(strcat('datafile_:_',datafiles(f).name,' is unreadable'));
                end
            end
            cd ../..
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %navigating to proper data directory
            runFiles = strcat('BB_',num2str(run));
            cd(runFiles);
            cd('fig_data');
            datafiles = dir(fullfile('*.txt'));
            %reading in data and storing to data matrix
            for k = 1:length(files2plot)
                dat = dlmread(files2plot(k)); 
                disp(strcat('Reading in : ',datafiles(k).name));
                data(a_a,a_b,run,k,:) = dat';
            end 
            cd ../..
        end
        cd ..
    end  
end
cd ../..
%end loading data

%replacing all zeros that appear after all cancer is killed with NaN values
% (so 0 values are not averaged with data from runs that kept going)
for i = 1:length(aa)
    for r = 1:bb
        for plt = 1:length(files2plot)
            lastnonzero = find(data(i,:,r,plt,:),1,'last');
            data(i,:,r,plt,lastnonzero:itters) = NaN;
        end
    end
end

%averaging data over runs and caculating CI
meanData = mean(data,3);
stdError = std(data,0,3)/sqrt(bb);
tScore = tinv([0.025  0.975],bb-1);
yCI95 = bsxfun(@times, stdError, tScore(1));

if runStats
    
    %ANOVA
    %running for stat sig on cell count 
    %for f = 1:length(data(1,1,1,:,1))
        
        %title(extractBefore(files2plot(f),'_'), 'fontsize', 18)
        
        %here 6 is the 6th datafile - numberofkills
        disp(extractBefore(files2plot(6),'_'));
        anovaData = squeeze(data(:,:,:,6,ANOVAsampleTime));
        [p, tbh, stats] = anova1(anovaData');
        [c,m,h,nms] = multcompare(stats);
    %end
    
    %Coralation
    
    
end

%Plotting data (averages)
%this plots one fig for each desired attribute, with a curve for each
%point in the param space, averaged over the runs.
if plotAvg
    %plotting param space averaged over runs
    %hold on
    for p = 1:length(data(1,1,1,:,1))
        figure
        title(extractBefore(files2plot(p),'_'), 'fontsize', 18)
        xlabel('Number of Iterations','fontsize', 14)
        ylabel(labels{p},'fontsize',14)
        %ylabel()
        hold on
        for AA = 1:length(aa)
            for AB = 1:length(ab)
                plot(squeeze(meanData(AA,AB,:,p,:)),'color',colors(AA,:),'LineWidth',2);
                
                
                if showErrorBars
                    errorbar(squeeze(meanData(AA,AB,:,p,:)),squeeze(yCI95(AA,AB,:,p,:)),'color',colors(AA,:),'HandleVisibility','off');
                end
                
            end
        end
        %generating figure labels
        led = [];
        for a = aa
            for b = ab
                led = [led,convertCharsToStrings(strcat(param1,'-',num2str(a),'--',param2,' - ',num2str(b)))];
            end
        end
        
        if contains(files2plot(p), 'Stems')
            legend(led,'Location','southwest')
        else
            legend(led,'Location','northwest')
        end
        if saveImgs
            cd Data
            mkdir SavedFigs
            cd SavedFigs
            mkdir(exType)
            cd(exType)
            
            print(strcat(extractBefore(files2plot(p),'.txt'),'.tiff'),'-dtiff',['-r' num2str(400)]);
            cd ../../..
        end
        hold off
        
    end
end

%plotting data (runs)
%this plots a fig for each point in the param space with a curve for each run 
if plotRuns
    for p = 1:length(data(1,1,1,:,1)) %figures
        for AA = 1:length(aa)
            for AB = 1:length(ab)
                figure
                xlabel('Number of Iterations','fontsize', 22)
                ylabel(replace(files2plot{p},'_','-'), 'fontsize', 12)
                hold on
                for run = 1:bb
                    plot(squeeze(data(AA,AB,run,p,:)),colors(AA,:),'DisplayName',strcat('AA_',AA,'__AB_','AB'));
                end
            end
        end
        %generating figure labels
        led = [];
        for a = aa
            for b = ab
                led = [led,convertCharsToStrings(strcat('AA-',num2str(a),'--','AB-',num2str(b)))];
            end
        end
        legend(led,'Location','southwest')
        if saveImgs
            cd Data
            mkdir SavedFigs
            cd SavedFigs
            mkdir(exType)
            cd(exType)
            saveas(gcf, strcat(extractBefore(files2plot(p),'.txt'),'.png'));
            cd ../..
        end
        hold off
    end
end
    
    



