function plotsummaryparams(mouseIDs,varargin)
%plot avg parameters (Ti,Te,VT,minVT,freq,peakflow) for each mouse, across
%conditions

%% create input parser

p=inputParser;

%define default conditions
defaultSaveCond=true;
defaultSampleRate=1000;
defaultPlotConditions={};
defaultSavePath='C:\Users\Iris Bachmutsky\Desktop\Bachmutskyetal.2020\In Vivo Example Breathing Data';
defaultSurgConditionOrder={'Baseline','PBC mOR KO','KF mOR KO','PBC KF mOR KO'};
defaultFolderName='Undefined Data';
defaultOxConditions={'normoxia','hypercapnia'};
defaultDrugConditions={'saline','morphine'};

%parse optional inputs
addRequired(p,'mouseIDs',@iscellstr);
addOptional(p,'SampleRate',defaultSampleRate,@isnumeric);
addOptional(p,'SaveCond',defaultSaveCond,@islogical);
addOptional(p,'PlotConditions',defaultPlotConditions,@iscellstr);
addOptional(p,'SavePath',defaultSavePath,@ischar);
addOptional(p,'SurgConditions',defaultSurgConditionOrder,@iscellstr);
addOptional(p,'FolderName',defaultFolderName,@ischar);
addOptional(p,'OxConditions',defaultOxConditions,@iscellstr);
addOptional(p,'DrugConditions',defaultDrugConditions,@iscellstr);

%parse
parse(p,mouseIDs,varargin{:});

%reset defaults
savecond=p.Results.SaveCond;
savepath=p.Results.SavePath;
surgcondorder=p.Results.SurgConditions;
fs=p.Results.SampleRate;
plotconditions=p.Results.PlotConditions;
numplots=length(plotconditions);
nummice=length(mouseIDs);
foldername=p.Results.FolderName;
oxconds=p.Results.OxConditions;
drugconds=p.Results.DrugConditions;

clear defaultFolderName defaultPlotConditions defaultSavePath defaultSurgConditionOrder defaultSaveCond defaultSampleRate varargin

%% generate distint colors, shapes for different conditions

tracecolors={'#cc0066','#cc5500','#00cc55','#0077cc','#9400D3'};
tracecolors=hex2rgb(tracecolors);
markertype={'d','s','o','^','v','d','s','o','^','v'};
grey=[187 187 187]/255;


%% set basic graphics properties for figure

set(groot, ...
    'DefaultFigureColor', 'w', ...
    'DefaultAxesLineWidth', 0.5, ...
    'DefaultAxesXColor', 'k', ...
    'DefaultAxesYColor', 'k', ...
    'DefaultAxesFontUnits', 'points', ...
    'DefaultAxesFontSize', 8, ...
    'DefaultAxesFontName', 'Helvetica', ...
    'DefaultLineLineWidth', 1, ...
    'DefaultTextFontUnits', 'Points', ...
    'DefaultTextFontSize', 8, ...
    'DefaultTextFontName', 'Helvetica', ...
    'DefaultAxesBox', 'off', ...
    'DefaultAxesTickLength', [0.02 0.025]);

% set the tickdirs to go out - need this specific order
set(groot, 'DefaultAxesTickDir', 'out');
set(groot, 'DefaultAxesTickDirMode', 'manual');

%% generate folder to save data in

cd(savepath)
if ~exist(strcat(savepath,'/',foldername),'dir')
    mkdir(foldername)
end
cd(foldername)
figpath=pwd;

%% set conditions to plot

% initialize structure
allmicevals=struct();

% will plot a figure for each oxygen condition:plotcondition pair
for oxindex=1:length(oxconds) %which oxygen condition is it (normoxia, hypercapnia
    oxcond=oxconds{oxindex};
    
    %initialize figure objects
    hfig = gobjects(numplots,1);
    
    for mouseindex=1:nummice
        mouseID=string(mouseIDs(mouseindex));
        
        %get basic info for where the mouse data is located
        [genotype,masterpathtodata,~,surgconds,expconds]=getmousepath(mouseID);
        
        cd(masterpathtodata);
        
        %% load in the data
        
        submouseid=strsplit(mouseID,'_');
        filename=strcat(submouseid{2},'_',submouseid{3},'_summary.mat');
        load(filename); %this is MouseSummary
        filename=strcat(submouseid{2},'_',submouseid{3},'_ExpParams.mat');
        load(filename); %this is ExpParams
        
        %% get names of all available files
        
        fnames_sum=fieldnames(MouseSummary);
        fnames_exp=fieldnames(ExpParams);
        
        %% get data for each plot condition, and fill in allmicevals
        
        flag=0;
        for plotindex=1:numplots
            plotcondition=plotconditions{plotindex};
            
            %initialize allmouse structure, which stores all data that will
            %be plotted in figures, ie. normalized plotcondition  data per
            %mouse
            header = string(surgcondorder);
            allmicevals.(plotcondition).(mouseID)(1,:)=header;
            allmicevals.(plotcondition).(mouseID)(2,:)=nan(1,length(surgcondorder));
            
            for surgindex=1:length(surgconds)
                surgcond=string(surgconds{surgindex}); %beforecre, aftercre, etc
                expcond=string(expconds(lower(surgconds)==lower(surgcond)));
                
                plotcondvals=struct();
                for drugindex=1:length(drugconds)
                    drugcond=drugconds{drugindex};
                    
                    %get index for data in that condition, in MouseSummary
                    searchterm=strcat(strrep(expcond,' ','_'),'_',oxcond,'_',drugcond);
                    isField=getmyfieldindex(fnames_sum,searchterm);
                    
                    %get the data
                    if any(isField)
                        tempsummary= MouseSummary.(fnames_sum{isField});
                    else
                        tempsummary = [];
                        flag=1;
                    end
                    
                    %get index for data in that condition, in ExpParams
                    isField=getmyfieldindex(fnames_exp,searchterm);
                    
                    %get the data
                    if any(isField)
                        tempexpparams=ExpParams.(fnames_exp{isField});
                    else
                        tempexpparams=[];
                        flag=2;
                    end
                    
                    %tell me if data is missing
                    if isempty(tempsummary) || isempty(tempexpparams) %if still not found display that it is not present
                        disp(sprintf('missing %s %s summary values',mouseID,searchterm))
                        continue
                    end
                    
                    % get the relevant parameter for each plot condition
                    switch plotcondition
                        case 'inspdur'
                            plotcondvals.(drugcond)=nanmean(tempsummary.inspDur); %in ms
                            normalization='divide'; %this defines whether data is normalized as morphine/saline, or morphine-saline
                        case 'instfreq'
                            breathstarts=tempsummary.breathStartInd; %in seconds
                            instfreq=breathstarts(2:end)-breathstarts(1:end-1);
                            instfreq=ones(length(instfreq),1)./instfreq;
                            plotcondvals.(drugcond)=nanmean([0 instfreq']);
                            normalization='divide';
                        case 'minVT'
                            plotcondvals.(drugcond)=nanmean(tempsummary.minVt);
                            normalization='divide';
                        case 'VT'
                            plotcondvals.(drugcond)=nanmean(tempsummary.Vt);
                            normalization='divide';
                        case 'peakflow'
                            plotcondvals.(drugcond)=nanmean(tempsummary.inspPeak);
                            normalization='divide';
                        case 'exppeak'
                            plotcondvals.(drugcond)=nanmean(tempsummary.expPeak);
                            normalization='divide';
                        case 'expdur'
                            plotcondvals.(drugcond)=nanmean(tempsummary.expDur);
                            normalization='divide';
                        case 'expdur_nopause'
                            plotcondvals.(drugcond)=nanmean(tempexpparams.Expdur);
                            normalization='divide';
                        case 'exppause'
                            plotcondvals.(drugcond)=nanmean(tempexpparams.Pause);
                            normalization='subtract';
                    end
                end
                
                %if you dont have a value for both saline and morphine,
                %move on
                if ~(isfield(plotcondvals,'saline') & isfield(plotcondvals,'morphine'))
                    disp(sprintf('missing %s %s summary values',mouseID,searchterm))                       
                    continue
                end
                
                %normalize by dividing by or subtracting saline condition
                if strcmp(normalization,'divide')
                    plotcondvals.normalized=plotcondvals.morphine/plotcondvals.saline;
                elseif strcmp(normalization,'subtract')
                    plotcondvals.normalized=plotcondvals.morphine-plotcondvals.saline;
                end
                
                %enter data into allmicevals
                plotlocation=find(strcmpi(surgcondorder,expcond)==1);
                allmicevals.(plotcondition).(mouseID)(2,plotlocation)=plotcondvals.normalized;
                
            end
            
        end
    end
    
    
    %% plot all values
    
    for plotindex=1:numplots
        plotcondition=plotconditions{plotindex};
        
        %initialize figure
        hfig(plotindex)=figure;
        title(plotcondition)
        
        mouseids_fortitle=[];
        %plot single mouse measurements
        for mouseindex=1:nummice
            mouseID=string(mouseIDs(mouseindex));
            mouseids_fortitle=[mouseids_fortitle mouseID];
            
            %ensure you are working on correct figure
            set(0, 'CurrentFigure', hfig(plotindex))
            
            %plot line connecting all values from one mouse
            plotvar=allmicevals.(plotcondition).(mouseID);
            tempx=1:length(surgcondorder);
            tempx=tempx(~ismissing(plotvar(2,:))); % in case there was on rec but not another
            tempy=str2double(plotvar(2,tempx));
            plot(tempx,tempy,'Color',grey,'LineWidth',1);
            hold on;
            
            %put marker on avg value from each mouse
            for ind=1:length(tempx)
                markerchoice=markertype{tempx(ind)};
                markercolor=tracecolors(tempx(ind),:);
                plot(tempx(ind),tempy(ind),'Marker',markerchoice,'MarkerEdgeColor',markercolor,'MarkerFaceColor',markercolor,'MarkerSize',10);
                hold on;
            end
            
        end
        
        % plot summary measurements across mice
        
        %set figure to top
        uistack(hfig(plotindex),'top');
        
        %merge data from all mice for the plot condition into a simplified matris
        %in order to get avg, stderr
        plotcondvals_allmice=allmicevals.(plotcondition);
        nummice=numel(fieldnames(plotcondvals_allmice));
        subfieldnames=fieldnames(plotcondvals_allmice);
        header =[' ' string(surgcondorder)];
        clear plotvar
        plotvar(1,:)=header;
        for ind=1:nummice
            subfieldname=string(subfieldnames(ind));
            plotvar(ind+1,:)=[subfieldname plotcondvals_allmice.(subfieldname)(2,:)];
        end
        
        %calculate avg and stderror for each condition, plot
        for ind=1:length(surgcondorder) %go through every conditiion
            data=str2double(plotvar(2:end,(ind+1)));
            data=data(~isnan(data));
            err=std(data)/sqrt(length(data)); %stderr
            errorbar(ind,mean(data),err,'o','MarkerEdgeColor','black','MarkerFaceColor','black','Color','black','LineStyle','none','LineWidth',1.5);
        end
        
        %% label & format your figure
        
        switch plotcondition
            case 'inspdur'
                ylabel('Normalized Inspiratory Duration');
            case 'instfreq'
                ylabel('Normalized Frequency');
            case 'minVT'
                ylabel('Normalized Minute Volume')
            case 'VT'
                ylabel('Normalized Tidal Volume')
            case 'peakflow'
                ylabel('Normalized Peak Flow')
            case 'exppeak'
                ylabel('Normalized Expiratory Peak Flow')
            case 'expdur'
                ylabel('Normalized Expiratory Duration')
                %             case 'expdurbyvolume'
                %                 ylabel('Normalized Expiratory Duration')
                %             case 'expdurbyvolume'
                %                 ylabel('Normalized Expiratory Duration 2')
            case 'expdur_nopause'
                ylabel('Normalized Expiratory Duration (w/o pause)')
            case 'exppause'
                ylabel('Normalized Expiratory Pause')
        end
        
        %label your axes
        xticks(1:length(surgcondorder));
        xticklabels(surgcondorder);
        xtickangle(45);
        plot([0 5],[1 1],'LineStyle','--','Color',grey); hold on;
        xlim([0 (length(surgcondorder)+1)])
        
        %add all mice into title, in order to keep track of which mice
        %you ran analysis on, move plotcondition to suptitle
        title(strcat(oxcond,' ',strrep(strjoin(mouseids_fortitle),'_',' ')));
        suptitle(plotcondition);
        
        %% save & close figure
        
        if savecond
            tempstr=strcat(figpath,'\',oxcond,'_',plotcondition,'_summary','.fig');
            savefig(hfig(plotindex),tempstr)
            tempstr=strcat(figpath,'\',oxcond,'_',plotcondition,'_summary','.pdf');
            saveas(hfig(plotindex),tempstr);
        end
        close(hfig(plotindex))
        
    end
end

%reset graphics settings
set(groot,'DefaultFigureColormap','remove')

end

function [isField]=getmyfieldindex(fnames,searchterm)

isField=contains(fnames,searchterm,'IgnoreCase',true);
notisField=contains(fnames,strcat('_',searchterm),'IgnoreCase',true);
isField=isField & ~notisField;
isField=find(isField,true,'last'); %in case there were multiple versions, shouldnt be

end