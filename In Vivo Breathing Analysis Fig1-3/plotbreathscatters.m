function plotbreathscatters(mouseID,varargin)
%generate scatterplots of multiple breath values

%%
%create input parser

p=inputParser;
%define default conditions
defaultSaveCond=true;
defaultSampleRate=1000;
defaultXCondition='Tidal Volume';
defaultYCondition='Instantaneous Frequency';
defaultdrugconds={'saline','morphine'};
defaultOxConditions={'normoxia','hypercapnia'};
defaultColor1={[0 0 0],[85,85,85]./255};
defaultColor2={[150,0,0]./255,[255,0,0]./255};

%parse optional inputs
addRequired(p,'animalID',@ischar);
addOptional(p,'SampleRate',defaultSampleRate,@isnumeric);
addOptional(p,'SaveCond',defaultSaveCond,@islogical);
addOptional(p,'XAxis',defaultXCondition,@ischar);
addOptional(p,'YAxis',defaultYCondition,@ischar);
addOptional(p,'DrugConditions',defaultdrugconds,@iscellstr);
addOptional(p,'OxConditions',defaultOxConditions,@iscellstr);
addOptional(p,'Color1',defaultColor1,@iscell);
addOptional(p,'Color2',defaultColor2,@iscell);

%parse
parse(p,mouseID,varargin{:});

%reset defaults
savecond=p.Results.SaveCond;
fs=p.Results.SampleRate;
plotonx=p.Results.XAxis;
plotony=p.Results.YAxis;
drugconds=p.Results.DrugConditions;
oxconds=p.Results.OxConditions;
color1=p.Results.Color1;
color2=p.Results.Color2;

%%
%get to the mouse data
[genotype,masterpathtodata,figpath,surgconds,expconds]=getmousepath(mouseID);

%%
%set figure folder

cd(figpath)
if ~exist(strcat(figpath,'/','scatterplots'),'dir')
    mkdir('scatterplots')
end
cd('scatterplots');
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
surgindex=1;
surgcond=string(surgconds{surgindex}); %inject1, inject2, etc
expcond=string(expconds(lower(surgconds)==lower(surgcond))); %ex: PBC MOR KO, PBC KF MOR KO, etc.

for oxindex=1:length(oxconds) %which breathing condition is it (normoxia, hypercapnia
    oxcond=oxconds{oxindex};
    
    %% initialize figure
    scatterfig=figure;
    
    ax1=axes;
    ax2=axes;
    
    for drugindex=1:length(drugconds)
        drugcond=drugconds{drugindex};
        plotcondition=strcat(strrep(mouseID,'_',' ')," ",expcond," ",oxcond);
        
        %%iterate through, pulling out the appropriate data
        searchterm=strcat(strrep(expcond,' ','_'),'_',oxcond,'_',drugcond);
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
        
        if strcmp(drugcond,'saline')
            
            %set color
            c1=color1{1};
            c2=color1{2};
            
            set(ax1,'Position',[.2 .15 .65 .815]);
            [grad,~]=colorGradient(c1,c2,length(xvar));
            scatter(ax1,xvar,yvar,sz,grad,'filled');
            hold on;
            colormap(ax1,grad);
            cb1 = colorbar(ax1,'Position',[0.04 .11 .0675 .815]); hold on;              cb1.TickLabels={'T0','Tn'};
            set(get(cb1,'title'),'string','saline','FontSize',10,'FontWeight','bold');
            
        elseif strcmp(drugcond,'morphine')
            
            %set color
            c1=color2{1};
            c2=color2{2};
            
            set(ax2,'Position',[.2 .15 .65 .815]);
            [grad,~]=colorGradient(c1,c2,length(xvar));
            scatter(ax2,xvar,yvar,sz,grad,'filled');
            hold on;
            colormap(ax2,grad);
            cb2 =colorbar(ax2,'Position',[.88 .11 .0675 .815]); hold on;
            cb2.Ticks=[0 1];
            cb2.TickLabels={'T0','Tn'};
            set(get(cb2,'title'),'string','morphine','FontSize',10,'FontWeight','bold');
        else
            error('you dun fucked up your drug conditions')
        end
        
    end
    
    %% set basic figure parameters
    
    set(0,'CurrentFigure',scatterfig);
    linkaxes([ax1,ax2])
    
    ax2.Visible = 'off';
    ax2.XTick = [];
    ax2.YTick = [];
    set(gca,'box','off')
    xlabel(ax1,xlab,'FontSize',10);
    ylabel(ax1,ylab,'FontSize',10);
    titlestr=oxcond;
    title(titlestr);
    suptitlestr=strcat(plotonx,' Vs. ',plotony);
    suptitle(suptitlestr);
    
    %% save figure
    if savecond
        tempstr=strcat(figpath,'\',mouseID,'_',oxcond,'_',strrep(suptitlestr," ",'_'),'.fig');
        savefig(scatterfig,tempstr)
        tempstr=strcat(figpath,'\',mouseID,'_',oxcond,'_',strrep(suptitlestr," ",'_'),'.pdf');
        saveas(scatterfig,tempstr);
    end
    close(scatterfig)
    
end
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
        plotlabel='Experiatory Pause (ms)';
end

end
