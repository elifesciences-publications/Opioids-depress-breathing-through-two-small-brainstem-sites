function plotavgnoverlaybreaths(mouseID,varargin)
%plot avg breath traces across conditions

%% create input parser

p=inputParser;

%define default conditions
defaultSaveCond=true;
defaultSampleRate=1000;
defaultSortCondition='inspiratory duration';
defaultLowerBound=20;
defaultUpperBound=40;
defaultPlotType='flow';
defaultOxConditions={'hypercapnia'};
defaultDrugConditions={'saline','morphine'};
defaultSavePath=nan;

%parse optional inputs
addRequired(p,'animalID',@ischar);
addOptional(p,'SampleRate',defaultSampleRate,@isnumeric);
addOptional(p,'SaveCond',defaultSaveCond,@islogical);
addOptional(p,'SortCondition',defaultSortCondition,@ischar);
addOptional(p,'LowerBound',defaultLowerBound,@isnumeric);
addOptional(p,'UpperBound',defaultUpperBound,@isnumeric);
addOptional(p,'PlotType',defaultPlotType,@ischar);
addOptional(p,'OxConditions',defaultOxConditions,@iscellstr);
addOptional(p,'DrugConditions',defaultDrugConditions,@iscellstr);
addOptional(p,'SavePath',defaultSavePath);

%parse
parse(p,mouseID,varargin{:});

%reset defaults
savecond=p.Results.SaveCond;
fs=p.Results.SampleRate;
sortcondition=p.Results.SortCondition;
lowerbound=p.Results.LowerBound;
upperbound=p.Results.UpperBound;
whichplot=p.Results.PlotType;
oxconds=p.Results.OxConditions;
drugconds=p.Results.DrugConditions;
savepath=p.Results.SavePath;

%% generate colors for plots
%for ex: all blues

tracecolors={'#99678D';'#906B92';'#866E96';'#7B7199';'#70749B';'#64779B';'#577A9A';'#4B7D99';'#3E7F95';'#338091';'#28828C';'#208386';'#1D847F';'#1F8478';'#258570';'#2D8568';'#36845F';'#408457';'#49834F';'#528247'};
tracecolors=hex2rgb(tracecolors);

%% get to the mouse data

[genotype,masterpathtodata,figpath,surgconds,expconds]=getmousepath(mouseID);

%% generate folder for data

cd(figpath)
if ~exist(strcat(figpath,'/','breath overlays'),'dir')
    mkdir('breath overlays')
end
cd('breath overlays');
figpath=pwd;
if ~exist(strcat(figpath,'/',sortcondition),'dir')
    mkdir(sortcondition)
end
cd(sortcondition)
figpath=pwd;

cd(masterpathtodata);

%% load in the data

submouseid=strsplit(mouseID,'_');
filename=strcat(submouseid{2},'_',submouseid{3},'_summary.mat');
load(filename); %this is MouseSummary
filename=strcat(submouseid{2},'_',submouseid{3},'_raw.mat');
load(filename); %this is AllBreaths

%% get names of all available files

fnames_sum=fieldnames(MouseSummary);
numfields_sum=length(fnames_sum);

fnames_breath=fieldnames(AllBreaths);
numfields_breath=length(fnames_breath);

%% iterate through conditions to plot

for oxindex=1:length(oxconds) %which breathing condition is it (normoxia, hypercapnia
    oxcond=oxconds{oxindex};
    
    for drugindex=1:length(drugconds)
        drugcond=drugconds{drugindex};
        
        % initialize figure
        suptitlestr=strcat(strrep(mouseID,'_',' ')," ",oxcond," ",drugcond);
        
        %flow
        if strcmp(whichplot,'all') || strcmp(whichplot,'flow')
            flowfig=figure;
            %tight_subplot(Nh, Nw, gap, marg_h, marg_w)
            flowaxes = tight_subplot(2,3,[.05 .05],[.15 .15],[.15 .15]);
            set(flowaxes,'box','off')
        end
        
        %volume
        if strcmp(whichplot,'all') || strcmp(whichplot,'volume')
            volumefig=figure;
            volumeaxes = tight_subplot(2,3,[.01 .03],[.15 .15],[.15 .15]);
            set(volumeaxes,'box','off')
        end
        
        %flow volume loop
        if strcmp(whichplot,'all') || strcmp(whichplot,'loop')
            loopfig=figure;
            loopaxes = tight_subplot(2,3,[.01 .03],[.15 .15],[.15 .15]);
            set(loopaxes,'box','off')
        end
        
        %% initialize variables for plotting
        missing=0;
        plotted=0;
        iter=1;
        fullfigmaxlength=0;
        fullfigmaxdepth_flow=0;
        fullfigmaxdepth_volume=0;
        
        for surgindex=1:length(surgconds)
            surgcond=string(surgconds{surgindex}); %beforecre, aftercre, etc
            expcond=string(expconds(lower(surgconds)==lower(surgcond)));
            fullcondition=strcat(strrep(mouseID,'_',' ')," ",expcond," ",oxcond," ",drugcond);
            
            %iterate through, pulling out the appropriate data
                        searchterm=strcat(strrep(expcond,' ','_'),'_',oxcond,'_',drugcond);
            [isField]=getmyfieldindex(fnames_sum,searchterm);
            
            %get the data
            if any(isField)
                tempsummary= MouseSummary.(fnames_sum{isField});
            else
                tempsummary = [];
            end
            
            %tell me if data is missing
            if isempty(tempsummary) %if still not found display that it is not present
                disp(sprintf('missing %s summary values',searchterm))
                missing=missing+1;
                continue
            end
            
            %get index for breath data in that condition, in AllBreahs
            [isField]=getmyfieldindex(fnames_breath,searchterm);

            %get the data
            if any(isField)
                tempbreaths= AllBreaths.(fnames_breath{isField});
            else
                tempbreaths = [];
            end
            
            %tell me if data is missing
            if isempty(tempbreaths) %if still not found display that it is not present
                disp(sprintf('missing %s breath values',searchterm))
                missing=missing+1;
                continue
            end
            
            if isempty(tempsummary) & isempty(tempbreaths)
                missing=missing-1;
            end
            
            
            %% define indices of breaths that fall into sort condition
    
            switch sortcondition
                case 'inspiratory duration'
                    var=tempsummary.inspDur;
                    metric='ms';
                case 'peak flow'
                    var=tempsummary.inspPeak;
                    metric='mls';
                case 'exp flow'
                    var=tempsummary.expPeak;
                    metric='mls';
                case 'tidal volume'
                    var=tempsummary.Vt;
                    metric='ml';
            end
            
            indices=find(var>=lowerbound & var<=upperbound);
            
            if isempty(indices)
                missing=missing+1;
                continue
            end
            
            %get approprate breaths (in flow, ml/sec)
            rawbreaths=tempbreaths(indices);
            
            %integrate for volume (ml)
            window=10; %10 ms window
            intbreaths=cellfun(@(data) (cumsum(data)./1000),rawbreaths,'UniformOutput',false);
            %             plotbreaths=plotbreaths(~cellfun(@isempty,plotbreaths));
            
            
            %% plot flow and volume, example traces of breaths
            
            maxplot=15;
            if length(rawbreaths)<maxplot maxplot=length(rawbreaths); end
            
            %define axes on plots
            maxlength=cellfun(@(x) length(x),rawbreaths);
            maxlength=max(maxlength); maxlength=maxlength(1);
            maxdepth_flow=cellfun(@(x) max(abs(x)),rawbreaths);
            maxdepth_flow=max(maxdepth_flow); maxdepth_flow=maxdepth_flow(1);
            maxdepth_volume=cellfun(@(x) max(abs(x)),intbreaths);
            maxdepth_volume=rmoutliers(maxdepth_volume); maxdepth_volume=max(maxdepth_volume); maxdepth_volume=maxdepth_volume(1);
            
            %% flow
            
            if strcmp(whichplot,'all') || strcmp(whichplot,'flow')
                %set which figure, for now
                set(0,'CurrentFigure',flowfig);
                
                %set plot axes
                axes(flowaxes(iter+3))
                
                %plot example traces
                for index=1:maxplot
                    temptime=1:1:length(rawbreaths{index}); %in ms
                    breath=rawbreaths{index}';
                    plot(temptime,breath,'k')
                    hold on;
                end
                set(gca,'box','off')
                xlim([0 maxlength]);
                ylim([-maxdepth_flow maxdepth_flow]);
                xlabel('ms')
                yticks([-20 -15 -10 -5 0 5 10 15 20]);
                
                %label yaxis for only the far left plot
                if iter==1
                    ylabel('Air Flow (ml/s)')
                end
            end
            %% volume
            
            if strcmp(whichplot,'all') || strcmp(whichplot,'volume')
                %set which figure, for now
                set(0,'CurrentFigure',volumefig);
                
                %set plot axes
                axes(volumeaxes(iter+3))
                
                %plot example traces
                for index=1:maxplot
                    breath=intbreaths{index}';
                    temptime=1:1:length(breath); %in ms
                    plot(temptime,breath,'k')
                    hold on;
                end
                set(gca,'box','off')
                xlim([0 maxlength]);
                ylim([-maxdepth_volume maxdepth_volume]);
                xlabel('ms')
                yticks(-1:0.25:1);
                
                %label yaxis for only the far left plot
                if iter==1
                    ylabel(' Volume (ml)')
                end
            end
            
            %% flow volume loop
            
            if strcmp(whichplot,'all') || strcmp(whichplot,'loop')
                %set which figure
                set(0,'CurrentFigure',loopfig);
                
                %set plot axes
                axes(loopaxes(iter+3))
                
                %plot example traces
                for index=1:maxplot
                    flowbreath=rawbreaths{index}';
                    volbreath=intbreaths{index}';
                    plot(volbreath,flowbreath,'k')
                    hold on;
                end
                set(gca,'box','off')
                xlim([-maxdepth_volume maxdepth_volume]);
                ylim([-maxdepth_flow maxdepth_flow]);
                xlabel('Volume (ml)')
                yticks([-20 -15 -10 -5 0 5 10 15 20]);
                
                %label yaxis for only the far left plot
                if iter==1
                    ylabel('Air Flow (ml/s)')
                end
            end
            
            plotted=plotted+1;
            
            %% plot flow and volume, avg traces of breaths
            
            %pad all breath waveforms to same length
            padrawbreaths=NaN(length(rawbreaths),maxlength);
            for index=1:length(rawbreaths)
                temptrace=rawbreaths{index};
                temptrace(end:maxlength)=0;
                padrawbreaths(index,:)=temptrace;
            end
            padintbreaths=NaN(length(intbreaths),maxlength);
            for index=1:length(intbreaths)
                temptrace=intbreaths{index};
                temptrace(end:maxlength)=0;
                padintbreaths(index,:)=temptrace;
            end
            
            %calculate average breaths
            meanrawbreaths=mean(padrawbreaths,1);
            stdrawbreaths=std(padrawbreaths,1);
            meanintbreaths=mean(padintbreaths,1);
            stdintbreaths=std(padintbreaths,1);
            
            %% flow
            
            if strcmp(whichplot,'all') || strcmp(whichplot,'flow')
                %set which figure
                set(0,'CurrentFigure',flowfig);
                
                %set plot axes
                axes(flowaxes(iter))
                
                %plot avg
                fillcolor={'#fc0303'};%{'#D3492B'};
                fillcolor=hex2rgb(fillcolor);
                temptime=1:1:maxlength;
                plot(temptime,meanrawbreaths,'k'); hold on;
                jbfill(temptime,meanrawbreaths,(meanrawbreaths+stdrawbreaths),fillcolor,[1 1 1],1,.7); hold on;
                jbfill(temptime,meanrawbreaths,(meanrawbreaths-stdrawbreaths),fillcolor,[1 1 1],1,.7); hold on;
                plot(temptime,meanrawbreaths,'k');
                
                %set axes params
                set(gca,'box','off')
                xlim([0 maxlength]);
                ylim([-maxdepth_flow maxdepth_flow]);
                set(gca,'XTickLabel',{[]});
                set(gca,'xtick',[])
                yticks([-20 -15 -10 -5 0 5 10 15 20]);
                
                %label yaxis for  only the far left plot
                if iter==1
                    ylabel('Air Flow (ml/s)')
                end
                title(expcond);
            end
            
            %% volume
            
            if strcmp(whichplot,'all') || strcmp(whichplot,'volume')
                %set which figure, for now
                set(0,'CurrentFigure',volumefig);
                
                %set plot axes
                axes(volumeaxes(iter))
                
                %plot avg
                fillcolor={'#5A24B3'};
                fillcolor=hex2rgb(fillcolor);
                temptime=1:1:maxlength;
                plot(temptime,meanintbreaths,'k'); hold on;
                jbfill(temptime,meanintbreaths,(meanintbreaths+stdintbreaths),fillcolor,[1 1 1],1,.7); hold on;
                jbfill(temptime,meanintbreaths,(meanintbreaths-stdintbreaths),fillcolor,[1 1 1],1,.7); hold on;
                plot(temptime,meanintbreaths,'k');
                
                %set axes params
                set(gca,'box','off')
                xlim([0 maxlength]);
                ylim([-maxdepth_volume maxdepth_volume]);
                set(gca,'XTickLabel',{[]});
                set(gca,'xtick',[])
                yticks(-1:0.25:1);
                %label yaxis for  only the far left plot
                if iter==1
                    ylabel('Volume (ml)')
                end
                title(expcond);
            end
            
            %% loop
            
            if strcmp(whichplot,'all') || strcmp(whichplot,'loop')
                %set which figure
                set(0,'CurrentFigure',loopfig);
                
                %set plot axes
                axes(loopaxes(iter))
               
                %%
                %plot avg
                fillcolor={'#002aff'};
                fillcolor=hex2rgb(fillcolor);
                plot(meanintbreaths,meanrawbreaths,'k'); hold on;
                %fill in std yaxis
                jbfill(meanintbreaths,meanrawbreaths,(meanrawbreaths+stdrawbreaths),fillcolor,[1 1 1],1,.7); hold on;
                jbfill(meanintbreaths,meanrawbreaths,(meanrawbreaths-stdrawbreaths),fillcolor,[1 1 1],1,.7); hold on;
                %fill in std xaxis
                % unsure how to do this
                plot(meanintbreaths,meanrawbreaths,'k'); hold on;
                
                %set axes params
                set(gca,'box','off')
                xlim([-maxdepth_volume maxdepth_volume]);
                ylim([-maxdepth_flow maxdepth_flow]);
                set(gca,'XTickLabel',{[]});
                set(gca,'xtick',[])
                yticks([-20 -15 -10 -5 0 5 10 15 20]);
                
                %label yaxis for  only the far left plot
                if iter==1
                    ylabel('Air Flow (ml/s)')
                end
                title(expcond);
            end
            
            %%
            
            %get maximal values for axes to reset all axes later
            if fullfigmaxlength<maxlength
                fullfigmaxlength=maxlength;
            end
            
            if fullfigmaxdepth_flow<maxdepth_flow
                fullfigmaxdepth_flow=maxdepth_flow;
            end
            
            if fullfigmaxdepth_volume<maxdepth_volume
                fullfigmaxdepth_volume=maxdepth_volume;
            end
            
            %next recording
            iter=iter+1;
            
        end
        
        %% reset axes parameters, save figures
        
        % decide if figures are worth saving (or if there is too much data
        %missing)
        if plotted<2
            close all
            continue;
        end
        if missing>2
            close all
            continue
        end
        
        %% flow
        
        if strcmp(whichplot,'all') || strcmp(whichplot,'flow')
            % reset xlim and ylim so they are the same bw all axes
            allAxesInFigure = findall(flowfig,'type','axes');
            for index=1:length(allAxesInFigure)
                ha=allAxesInFigure(index);
                set(ha,'xlim',[0 fullfigmaxlength])
                set(ha,'ylim',[-fullfigmaxdepth_flow fullfigmaxdepth_flow])
            end
            % assign figure title
            set(0,'CurrentFigure',flowfig);
            suptitle(strcat(suptitlestr," ",num2str(lowerbound),metric,'-',num2str(upperbound),metric));
        end
        
        %% volume
        
        if strcmp(whichplot,'all') || strcmp(whichplot,'volume')
            % reset xlim and ylim so they are the same bw all axes
            allAxesInFigure = findall(volumefig,'type','axes');
            for index=1:length(allAxesInFigure)
                ha=allAxesInFigure(index);
                set(ha,'xlim',[0 fullfigmaxlength])
                set(ha,'ylim',[-fullfigmaxdepth_volume fullfigmaxdepth_volume])
            end
            % assign figure title
            set(0,'CurrentFigure',volumefig);
            suptitle(strcat(suptitlestr," ",num2str(lowerbound),metric,'-',num2str(upperbound),metric));
        end
        
        %% loop
                
        if strcmp(whichplot,'all') || strcmp(whichplot,'loop')
            % reset xlim and ylim so they are the same bw all axes
            allAxesInFigure = findall(loopfig,'type','axes');
            for index=1:length(allAxesInFigure)
                ha=allAxesInFigure(index);
                set(ha,'xlim',[-maxdepth_volume maxdepth_volume])
                set(ha,'ylim',[-maxdepth_flow maxdepth_flow])
            end
            % assign figure title
            set(0,'CurrentFigure',loopfig);
            suptitle(strcat(suptitlestr," ",num2str(lowerbound),metric,'-',num2str(upperbound),metric));
        end      
        
        %% save plots
        
        if savecond
            
            cd(figpath)
            if strcmp(whichplot,'all') || strcmp(whichplot,'flow')
                savefig(flowfig,strcat(figpath,'\',suptitlestr,'_',sortcondition,'_',num2str(lowerbound),metric,'_',num2str(upperbound),metric,'_flow.fig'))
                saveas(flowfig,strcat(figpath,'\',suptitlestr,'_',sortcondition,'_',num2str(lowerbound),metric,'_',num2str(upperbound),metric,'_flow.pdf'));
                close(flowfig)
            end
            if strcmp(whichplot,'all') || strcmp(whichplot,'volume')
                savefig(volumefig,strcat(figpath,'\',suptitlestr,'_',sortcondition,'_',num2str(lowerbound),metric,'_',num2str(upperbound),metric,'_vol.fig'))
                saveas(volumefig,strcat(figpath,'\',suptitlestr,'_',sortcondition,'_',num2str(lowerbound),metric,'_',num2str(upperbound),metric,'_vol.pdf'))
                close(volumefig)
            end
            if strcmp(whichplot,'all') || strcmp(whichplot,'loop')
                savefig(loopfig,strcat(figpath,'\',suptitlestr,'_',sortcondition,'_',num2str(lowerbound),metric,'_',num2str(upperbound),metric,'_loop.fig'))
                saveas(loopfig,strcat(figpath,'\',suptitlestr,'_',sortcondition,'_',num2str(lowerbound),metric,'_',num2str(upperbound),metric,'_loop.pdf'))
                close(loopfig)
            end
        end
        
    end
end

end

function [isField]=getmyfieldindex(fnames,searchterm)

isField=contains(fnames,searchterm);
notisField=contains(fnames,strcat('_',searchterm));
isField=isField & ~notisField;
isField=find(isField,true,'last'); %in case there were multiple versions, shouldnt be

end
