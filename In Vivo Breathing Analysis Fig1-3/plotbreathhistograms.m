function plotbreathhistograms(mouseID,varargin)
%generate histograms of breath parameters across different recording
%conditions


%% create input parser

p=inputParser;

%define default conditions
defaultSaveCond=true;
defaultSampleRate=1000;
defaultPlotCondition='inspiratory duration';
defaultOxConditions={'normoxia','hypercapnia'};
defaultDrugConditions={'morphine'};

%parse optional inputs
addRequired(p,'animalID',@ischar);
addOptional(p,'SaveCond',defaultSaveCond,@islogical);
addOptional(p,'SampleRate',defaultSampleRate,@isnumeric);
addOptional(p,'PlotCondition',defaultPlotCondition,@ischar);
addOptional(p,'OxConditions',defaultOxConditions,@iscellstr)
addOptional(p,'DrugConditions',defaultDrugConditions,@iscellstr);

%parse
parse(p,mouseID,varargin{:});

%set vars
savecond=p.Results.SaveCond;
fs=p.Results.SampleRate;
plotcondition=p.Results.PlotCondition;
oxconds=p.Results.OxConditions;
drugconds=p.Results.DrugConditions;


%% generate colors for plots (array of distinct colors to plot)

tracecolors={'#812ED1','#45FF00','#D538FF','#00BDFF','#FF5F00','#143CCC','#FFD801'};
tracecolors=hex2rgb(tracecolors);

%% set up figure folder for each condition

[genotype,masterpathtodata,figpath,surgconds,expconds]=getmousepath(mouseID);

cd(figpath)
if ~exist(strcat(figpath,'/','histograms'),'dir')
    mkdir('histograms')
end
cd('histograms')
figpath=pwd;

if ~exist(strcat(figpath,'/',plotcondition),'dir')
    mkdir(plotcondition)
end
cd(plotcondition)
figpath=pwd;

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

%% iterate through, pulling out the appropriate data

for oxindex=1:length(oxconds) %which breathing condition is it (normoxia, hypercapnia
    oxcond=oxconds{oxindex};
    
    for drugindex=1:length(drugconds)
        drugcond=drugconds{drugindex};%which druc condition (saline, morphine etc)
        
        %initialize figure
        histfig=figure;
        varSummary=[];
        legendsum=[];
        
        for surgindex=1:length(surgconds)
            surgcond=string(surgconds{surgindex}); %which injection (ex: injection 1)
            expcond=string(expconds(lower(surgconds)==lower(surgcond))); %what experimental condition (ex: PBC MOR KO)
            
            %% get the data you want, from the two summary data matrices 
            
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
            
            %get data based on condition set for plotting
            var=getmydata(tempsummary,tempexpparams,plotcondition);
            if isempty(var)
                flag=1;
                continue
            end
            
            varSummary{end+1}=var;
            legendsum{end+1}=expcond;
            
            
        end
        
        %if this condition was never recorded
        if isempty(varSummary)
            flag=1;
            continue
        end
        
        %% add a control trace for comparison on morphine histograms only
      
        %either combine all saline conditions together
           combinesaline=true;
        % or seperately plot each one
%         combinesaline=false;
        if strcmp(drugcond,'morphine')  || strcmp(drugcond,'fentanyl')
            
            controlvarSummary=[];
            controllegendSum=[];
            for surgindex=1:length(surgconds)
                surgcond=string(surgconds{surgindex}); %beforecre, aftercre, etc
                expcond=string(expconds(lower(surgconds)==lower(surgcond)));
                
                %get index for data in that condition, in MouseSummary
                searchterm=strcat(strrep(expcond,' ','_'),'_',oxcond,'_','saline');
                
                isField=getmyfieldindex(fnames_sum,searchterm);
                
                if any(isField)
                    tempsummary= MouseSummary.(fnames_sum{isField});
                else
                    disp(sprintf('missing %s %s summary values',mouseID,searchterm))
                    tempsummary=[];
                    flag=2;
                    continue
                end
                
                isField=getmyfieldindex(fnames_exp,searchterm);
                
                if any(isField)
                    tempexpparams= ExpParams.(fnames_exp{isField});
                else
                    disp(sprintf('missing %s %s summary values',mouseID,searchterm))
                    tempsummary=[];
                    flag=2;
                    continue
                end
                
                %get data based on condition set for plotting
                var=getmydata(tempsummary,tempexpparams,plotcondition);
                if isempty(var)
                    continue
                end
                
                controlvarSummary{end+1}=var;
                controllegendSum{end+1}=strcat('Saline'," ",expcond);
            end
            
            if combinesaline==true
                try
                    combination=cell2mat(controlvarSummary')';
                catch
                    combination=cell2mat(controlvarSummary);
                end
                varSummary{end+1}=combination;
                legendsum{end+1}='All Saline';
            else
                varSummary=[varSummary controlvarSummary];
                legendsum=[legendsum controllegendSum];
            end
            
        end
        
        
        %% plot all data
        
        legendsum=cellfun(@char,legendsum,'UniformOutput',false);
        
        %plot all values
        colormap(tracecolors);
        %for exppause plot normal, and logorithmic scale histogram
        if strcmp(plotcondition,'exppause')
            subplot(1,2,1)
            nhist(varSummary,'samebins','pdf','smooth','legend',legendsum,'location','east'); hold on;
            xlim([-100 1500])
            xlabel('ms')
            box off
            subplot(1,2,2)
            logVarSummary=cellfun(@log10, varSummary,'UniformOutput',false);            
            nhist(logVarSummary,'maxbins',50,'samebins','pdf','smooth','legend',legendsum,'location','east'); hold on;
            xlim([0 log10(1500)])
            xlabel('log(expiratory pause)')
        else
            nhist(varSummary,'samebins','pdf','smooth','legend',legendsum,'location','east'); hold on;
            %label your figure
            switch plotcondition
                case 'inspiratory duration'
                    xlabel('ms');
                case 'instantaneous frequency'
                    xlabel('Hz');
                case 'expiratory duration'
                    xlabel('ms');
                case 'inspiratory peak'
                    xlabel('ml/s');
                case 'expiratory peak'
                    xlabel('ml/s');
                case 'tidal volume'
                    xlabel('ml');
                case 'expdur_nopause'
                    xlabel('ms');
                case 'exppause'
                    xlabel('ms');
            end
        end
        
        titlestr=strcat(strrep(mouseID,'_'," ")," ",oxcond," ",drugcond);
        suptitle(titlestr)
        
        if savecond
            
            set(0,'CurrentFigure',histfig);
            tempstr=strcat(figpath,'\',strrep(titlestr," ",'_'),'_',strrep(plotcondition,' ','_'),'.fig');
            savefig(histfig,tempstr)
            tempstr=strcat(figpath,'\',strrep(titlestr," ",'_'),'_',strrep(plotcondition,' ','_'),'.pdf');
            saveas(histfig,tempstr);
            
        end
        close(histfig)
        
    end
    
end
end

function [isField]=getmyfieldindex(fnames,searchterm)

isField=contains(fnames,searchterm);
notisField=contains(fnames,strcat('_',searchterm));
isField=isField & ~notisField;
isField=find(isField,true,'last'); %in case there were multiple versions, shouldnt be

end

function [var]=getmydata(summary,expparams,plotcondition)

var=[];
switch plotcondition
    case 'inspiratory duration'
        if ~isfield(summary,'inspDur')
            return
        end
        var=summary.inspDur;
    case 'instantaneous frequency'
        if ~isfield(summary,'breathStartInd')
            return
        end
        breathstarts=summary.breathStartInd; %in seconds
        instfreq=breathstarts(2:end)-breathstarts(1:end-1);
        instfreq=ones(length(instfreq),1)./instfreq; 
        var=[0 instfreq'];
    case 'expiratory duration'
        if ~isfield(summary,'expDur')
            return
        end
        var=summary.expDur(:,1);
    case 'inspiratory peak'
        if ~isfield(summary,'inspPeak')
            return
        end
        var=summary.inspPeak(:,1);
    case 'expiratory peak'
        if ~isfield(summary,'expPeak')
            return
        end
        var=summary.expPeak(:,1);
    case 'tidal volume'
        if ~isfield(summary,'Vt')
            return
        end
        var=summary.Vt;
    case 'expdur_nopause'
        if ~isfield(expparams,'Expdur')
            return
        end
        var=expparams.Expdur;
    case 'exppause'
        if ~isfield(expparams,'Pause')
            return
        end
        var=expparams.Pause;
end

end

