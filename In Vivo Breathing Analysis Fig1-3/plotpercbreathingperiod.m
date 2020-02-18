function plotpercbreathingperiod(mouseID,varargin)
%the purpose of this code is to plot the relative contribution of each
%breath period on the total breathing time

%%
%create input parser

p=inputParser;

%define default conditions
defaultSaveCond=true;
defaultSampleRate=1000;
defaultOxConditions={'hypercapnia','normoxia'};
defaultDrugConditions={'saline','morphine'};
defaultSavePath='C:\Users\Iris Bachmutsky\Desktop\TempFigs';


%parse optional inputs
addRequired(p,'mouseID',@ischar);
addOptional(p,'SaveCond',defaultSaveCond,@islogical);
addOptional(p,'SampleRate',defaultSampleRate,@isnumeric);
addOptional(p,'OxConditions',defaultOxConditions,@iscellstr);
addOptional(p,'DrugConditions',defaultDrugConditions,@iscellstr);
addOptional(p,'SavePath',defaultSavePath,@ischar);

%parse
parse(p,mouseID,varargin{:});

%reset defaults
savecond=p.Results.SaveCond;
fs=p.Results.SampleRate;
oxconds=p.Results.OxConditions;
drugconds=p.Results.DrugConditions;
savepath=p.Results.SavePath;

%% set figure and data paths

[genotype,masterpathtodata,figpath,surgconds,expconds]=getmousepath(mouseID);

%% generate folder for data

cd(figpath)
if ~exist(strcat(figpath,'/','pause'),'dir')
    mkdir('breath period')
end
cd('breath period');
figpath=pwd;

%% load in the data

cd(masterpathtodata);
submouseid=strsplit(mouseID,'_');
filename=strcat(submouseid{2},'_',submouseid{3},'_summary.mat');
load(filename); %this is MouseSummary

filename=strcat(submouseid{2},'_',submouseid{3},'_ExpParams.mat');
load(filename);

%% get names of all available files

fnames_sum=fieldnames(MouseSummary);
fnames_exp=fieldnames(ExpParams);

%% set conditions to plot

numcols=length(drugconds);

flag=0;
for surgindex=1:length(surgconds)
    surgcond=string(surgconds{surgindex}); %beforecre, aftercre, etc
    expcond=string(expconds(lower(surgconds)==lower(surgcond)));
    
    freqcomp=[];
    for oxcondindx=1:length(oxconds)
        oxcond=oxconds{oxcondindx};
        fullcondition=strcat(mouseID,'_',expcond,'_',oxcond);
        
        %initialized figures
        piefig=figure;
        set(gcf,'Position',[100 100 600 300])
        piefig.Units='normalized';
        barfig=figure;
        axes1=subplot(1,3,1);
        axes2=subplot(1,3,2);
        axes3=subplot(1,3,3);
        
        %initialize data arrays
        insparray=NaN(1,length(drugconds));
        exparray=NaN(1,length(drugconds));
        pausearray=NaN(1,length(drugconds));
        
        for drugcondindx=1:length(drugconds)
            drugcond=drugconds(drugcondindx);
            
            %get the data you want
            searchterm=strcat(strrep(expcond,' ','_'),'_',oxcond,'_',drugcond);
            [isField]=getmyfieldindex(fnames_sum,searchterm);
           
            if any(isField)
                tempsummary= MouseSummary.(fnames_sum{isField});
            else
                disp(sprintf('missing %s summary values',searchterm))
                flag=1;
            end
            
            [isField]=getmyfieldindex(fnames_exp,searchterm);
           
            if any(isField)
                tempexpparams= ExpParams.(fnames_exp{isField});
            else
                disp(sprintf('missing %s summary values',searchterm))
                flag=1;
            end
            
            %get inspdur and expdur for each breath
            inspdur=tempsummary.inspDur;
            expdur=tempsummary.expDur;
            
            %get ExpParams values for expiratory duration and pause fore
            %each breath
            expdur_v2=tempexpparams.Expdur;
            pause=tempexpparams.Pause;
            
            totalbreathdur=inspdur+expdur;
            
            if strcmp(drugcond,'saline')
                breathstarts=tempsummary.breathStartInd;
                breathstarts=breathstarts(~isnan(breathstarts));
                freqcomp.saline=length(breathstarts)/(breathstarts(end)-breathstarts(1));
            elseif strcmp(drugcond,'morphine')
                breathstarts=tempsummary.breathStartInd;
                breathstarts=breathstarts(~isnan(breathstarts));
                freqcomp.morphine=length(breathstarts)/(breathstarts(end)-breathstarts(1));
            end
            
            %ensure expiratory data is matching up to expectation
            if ~isequal(expdur,(expdur_v2+pause))
                disp(sprintf('problem with %s expiratory values',searchterm))
            end
            
            %% set up data for plots
            
            avbreathlength=sum(totalbreathdur)/length(totalbreathdur);
            avinspdur=sum(inspdur)/length(inspdur);
            avexpdur=sum(expdur_v2)/length(expdur_v2);
            avpausedur=sum(pause)/length(pause);
            
            index=find(strcmp(drugconds,drugcond));
            insparray(index)=avinspdur;
            exparray(index)=avexpdur;
            pausearray(index)=avpausedur;
            
            
            %% plot data for pie figure
            
            set(0,'CurrentFigure',piefig)
            %normalized to either recording duration
            %pievals=[ sum(inspdur)/sum(totalbreathdur) sum(expdur_v2)/sum(totalbreathdur) sum(pause_v2)/sum(totalbreathdur)];
            %or to num breaths
            pievals=[avinspdur/avbreathlength avexpdur/avbreathlength avpausedur/avbreathlength];
            percvals=100*pievals;
            labels={sprintf('Insp (%0.1f%%)',percvals(1)),sprintf('Exp (%0.1f%%)',percvals(2)),sprintf('Pause (%0.1f%%)',percvals(3))};
            numplot=drugcondindx;
            hold on;
            subplot(1,numcols,numplot);
            h=get(gca);
            set(gca,'Position',[0.15+(numplot-1)/2.5 0.3 1/3 0.5]);  % This I what I added, You need to play with this
            explode=[0 0 1];
            pie(gca,pievals,explode,labels);
            set(gca,'xcolor','none')
            set(gca,'ycolor','none')
            title(drugcond, 'FontSize',12,'Units', 'normalized', 'Position', [0.5, -0.3, 0]); % Set Title with correct Position
            
        end

        if(flag==1)
            flag=0;
            close all;
            break
        end

        %% plot bar figure
        
        set(0,'CurrentFigure',barfig)
        xarray=1:length(drugconds);
        grey=[85,85,85]./255;
        bar(axes1,xarray,insparray,'FaceColor',grey);
        bar(axes2,xarray,exparray,'FaceColor',grey);
        bar(axes3,xarray,pausearray,'FaceColor',grey);
        labels=drugconds;
        setbasicfigparams(axes1,labels);
        title(axes1,'Inspiration');
        setbasicfigparams(axes2,labels);
        title(axes2,'Expiration');
        setbasicfigparams(axes3,labels);
        title(axes3,'Pause');
        
        suptitle('Average Changes to Breath Periods');
        
        %% any extra formatting to piefig, save
        set(0,'CurrentFigure',piefig)
        if ~isempty(freqcomp)
            sgtitle(sprintf('%s Rate Depression Ratio: %0.1f',expcond,freqcomp.morphine/freqcomp.saline));
        end
        
        if savecond
            %             cd(savepath)
            cd(figpath)
            set(0,'CurrentFigure',piefig)
            savefig(piefig,strcat(fullcondition,'_Pie.fig'))
            saveas(piefig,strcat(fullcondition,'_Pie.pdf'));
            close(piefig)
            set(0,'CurrentFigure',barfig)
            savefig(barfig,strcat(fullcondition,'_Bar.fig'));
            saveas(barfig,strcat(fullcondition,'_Bar.pdf'));
            close(barfig)
        end
    end
end
end

function setbasicfigparams(ha,labels)
set(ha,'ylim',[0 500])
ylabel(ha,'ms');
xtickarray=1:length(labels);
xticks(ha,xtickarray);
xticklabels(ha,labels);
xtickangle(ha,45);
xlim(ha,[0 3]);
set(ha,'box','off');
end

function [isField]=getmyfieldindex(fnames,searchterm)

isField=contains(fnames,searchterm);
notisField=contains(fnames,strcat('_',searchterm));
isField=isField & ~notisField;
isField=find(isField,true,'last'); %in case there were multiple versions, shouldnt be

end


