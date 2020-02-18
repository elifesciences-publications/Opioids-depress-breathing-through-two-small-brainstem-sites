function plotbreathscatters(mouseID,varargin)
%generate scatterplots of multiple breath values, inputting which
%conditions you want to plot

%% create input parser

p=inputParser;
%define default conditions
defaultSaveCond=true;
defaultSampleRate=1000;
defaultXCondition='Tidal Volume';
defaultYCondition='Instantaneous Frequency';
defaultPlotConditions={'Baseline_normoxia_saline','Baseline_normoxia_morphine'};
defaultColorMatrix= {[0 0 0],[85,85,85]./255;[150,0,0]./255,[255,0,0]./255;};
defaultSubFilename='PBC mOR KO';

%parse optional inputs
addRequired(p,'animalID',@ischar);
addOptional(p,'SampleRate',defaultSampleRate,@isnumeric);
addOptional(p,'SaveCond',defaultSaveCond,@islogical);
addOptional(p,'XAxis',defaultXCondition,@ischar);
addOptional(p,'YAxis',defaultYCondition,@ischar);
addOptional(p,'PlotConditions',defaultPlotConditions,@iscellstr);
addOptional(p,'ColorMatrix',defaultColorMatrix,@iscell);
addOptional(p,'SubFilename',defaultSubFilename,@ischar);

%parse
parse(p,mouseID,varargin{:});

%reset defaults
savecond=p.Results.SaveCond;
fs=p.Results.SampleRate;
plotonx=p.Results.XAxis;
plotony=p.Results.YAxis;
PlotConditions=p.Results.PlotConditions;
ColorMatrix=p.Results.ColorMatrix;
subfilename=p.Results.SubFilename;

if length(PlotConditions)~=length(ColorMatrix)
    error('Mismatch between plotting variables and colors')
end

%%
%get to the mouse data
[genotype,masterpathtodata,figpath,surgconds,expconds]=getmousepath(mouseID);

%%
%set figure folder

cd(figpath)
if ~exist(strcat(figpath,'/','scatterplots'),'dir')
    mkdir('scatterplots')
end
cd('scatterplots')
figpath=pwd;

if ~exist(strcat(figpath,'/',subfilename),'dir')
    mkdir(subfilename)
end
cd(subfilename)
figpath=pwd;

cd(masterpathtodata);

%%
%load in the data

submouseid=strsplit(mouseID,'_');
filename=strcat(submouseid{2},'_',submouseid{3},'_summary.mat');
load(filename); %this is MouseSummary

filename=strcat(submouseid{2},'_',submouseid{3},'_ExpParams.mat');
load(filename); %this is ExpParams

%%
%get names of all available files

fnames_sum=fieldnames(MouseSummary);
fnames_exp=fieldnames(ExpParams);

%% iterate through conditions to plot

%initialize conditions
missing=0;
numconditions=length(PlotConditions);

%% initialize figure
scatterfig=figure;
axarray=[];
i=1;
while i<=numconditions
    newaxes=addaxes(axarray,scatterfig);
    axarray=newaxes;
    i=i+1;
end

for plotcondindex=1:length(PlotConditions)
    plotcondition=PlotConditions{plotcondindex};
    color=ColorMatrix(plotcondindex,:);
    currentaxes=axarray(plotcondindex);
    
    allconds=strsplit(plotcondition,'_');
    expcond=strjoin(allconds(1:end-1));
    oxcond=allconds{end-1};
    drugcond=allconds{end};
    clear allconds
    
    %%iterate through, pulling out the appropriate data
    searchterm=plotcondition;
    isField=getmyfieldindex(fnames_sum,searchterm);
    
    if any(isField)
        tempsummary= MouseSummary.(fnames_sum{isField});
    else
        disp(sprintf('missing %s %s summary values',mouseID,searchterm))
        tempsummary=[];
        flag=1;
        continue
    end
    
    isField=getmyfieldindex(fnames_exp,searchterm);
    
    if any(isField)
        tempexpparams= ExpParams.(fnames_exp{isField});
    else
        disp(sprintf('missing %s %s summary values',mouseID,searchterm))
        tempsummary=[];
        flag=1;
        continue
    end
    
    %tell me if data is missing
    if isempty(tempsummary) %if still not found display that it is not present
        disp(sprintf('missing %s summary values',searchterm))
        missing=missing+1;
        continue
    end
    
    %% pick which condition to plot on xaxis
    
    [xvar,xlab]= getmydata(tempsummary,tempexpparams,plotonx);
    
    %% pick which condition to plot on yaxis
    
    [yvar,ylab]= getmydata(tempsummary,tempexpparams,plotony);
    
    %% check values are appropriate size
    
    if isempty(xvar) || isempty(yvar)
        error('Incorrect Assignment of Values')
    elseif length(xvar)~=length(yvar)
        error('Axis Values are Mismatched')
    end
    
    %% decide on color palette, plot values
    
    set(0,'CurrentFigure',scatterfig);
    sz = 25;
    
    %set color
    c1=color{1};
    c2=color{2};
    
    uistack(currentaxes, 'top')
    set(currentaxes,'Position',[.2 .15 .65 .815]);
    [grad,~]=colorGradient(c1,c2,length(xvar));
    scatter(currentaxes,xvar,yvar,sz,grad,'filled');
    hold on;
    
end

%% set basic figure parameters

set(0,'CurrentFigure',scatterfig);
linkaxes(axarray)

i=2;
while i<=numconditions
    axarray(i).Visible = 'off';
    axarray(i).Visible='off';
    axarray(i).XTick=[];
    axarray(i).YTick=[];
    set(axarray(i),'box','off')
%     currentaxes=axarray(i);
%     currentaxes.Visible='off';
%     currentaxes.XTick=[];
%     currentaxes.YTick=[];
%     set(currentaxes,'box','off')
    i=i+1;
end

xlim([0 30]);
ylim([-20 0]);
xlabel(axarray(1),xlab,'FontSize',10);
ylabel(axarray(1),ylab,'FontSize',10);
titlestr=strcat(strrep(mouseID,'_',' ')," ",oxcond);
set(get(axarray(1),'title'),'string',titlestr);

suptitlestr=strcat(plotonx,' Vs. ',plotony);
suptitle(suptitlestr);

%% save figure
if savecond
    tempstr=strcat(figpath,'\',strrep(titlestr," ",'_'),'_',strrep(suptitlestr," ",'_'),'.fig');
    savefig(scatterfig,tempstr)
    tempstr=strcat(figpath,'\',strrep(titlestr," ",'_'),'_',strrep(suptitlestr," ",'_'),'.pdf');
    saveas(scatterfig,tempstr);
end
close(scatterfig)

end


function [isField]=getmyfieldindex(fnames,searchterm)

isField=contains(fnames,searchterm);
notisField=contains(fnames,strcat('_',searchterm));
isField=isField & ~notisField;
isField=find(isField,true,'last'); %in case there were multiple versions, shouldnt be

end

function [plotvar,plotlabel]= getmydata(summary,expsummary,plotcondition)

plotvar=[];
plotlabel=[];
switch plotcondition
    case 'Inspiratory Duration'
        if ~isfield(summary,'inspDur') %fixed this in gen code
            return
        end
        plotvar=summary.inspDur; %in ms
        plotlabel='Inspiratory Duration (ms)';
    case 'Instantaneous Frequency'
        if ~isfield(summary,'breathStartInd')
            return
        end
        breathstarts=summary.breathStartInd; %in seconds
        instfreq=breathstarts(2:end)-breathstarts(1:end-1);
        instfreq=ones(length(instfreq),1)./instfreq;
        plotvar=instfreq';
        plotlabel='Instantaneous Frequency (Hz)';
    case 'Expiratory Duration'
        if ~isfield(summary,'expDur') %fixed this in gen code
            return
        end
        plot var=summary.expDur(:,1);
        plotlabel='Expiratory Duration (ms)';
    case 'Peak Flow'
        if ~isfield(summary,'inspPeak') %fixed this in gen code
            return
        end
        plotvar=summary.inspPeak(:,1);
        plotlabel='Peak Inspiratory Flow (ml/s)';
    case 'Peak Exp Flow'
        if ~isfield(summary,'expPeak') %fixed this in gen code
            return
        end
        plotvar=summary.expPeak(:,1);
        plotlabel='Peak Expiratory Flow (ml/s)';
    case 'Tidal Volume'
        if ~isfield(summary,'Vt') %fixed this in gen code
            return
        end
        plotvar=summary.Vt;
        plotlabel='Tidal Volume (ml)';
    case 'Expiratory Duration w/o Pause'
        if ~isfield(expsummary,'Expdur_v2')
            return
        end
        plotvar=expsummary.Expdur_v2;
        plotlabel='Expiratory Duration (ms)';
    case 'Expiratory Pause'
        if ~isfield(expsummary,'Pause_v2')
            return
        end
        plotvar=expsummary.Pause_v2;
        plotlabel='Expiratory Pause (ms)';
    case 'Net Volume'
         if ~isfield(summary,'breathVt')
             return
         end
         plotvar=summary.breathVt;
         plotlabel='Net Volume (ml)'
    case 'Expiratory Volume'
        if ~isfield(summary,'expVt')
            return
        end
        plotvar=summary.expVt;
        plotlabel='Expiratory Volume (ml)';
    case 'Expiratory Volume (no pause)'
        if ~isfield(expsummary,'ExpVol')
            return
        end
        plotvar=expsummary.ExpVol;
        plotlabel='Expiratory Volume (no pause) (ml)';
    case 'Pause Volume'
        if ~isfield(expsummary,'PauseVol')
            return
        end
        plotvar=expsummary.PauseVol;
        plotlabel='Pause Volume (ml)';
end

end

function newaxarray=addaxes(axarray,figname)

if isempty(axarray)
    newaxarray=axes(figname);
else
    newaxes=axes(figname);
    newaxarray=[axarray newaxes];
end

end
