function plotpercalive(DoseResponseCurve,genotypes,dose,varargin)
% plot bar graph of which percent of slices are dead at relative dose per
% genotype

%%
%create input parser
p=inputParser;

%define default conditions
defaultFigurePath='C:\Users\Iris Bachmutsky\Desktop\Bachmutskyetal.2020\In Vitro Example Extracellular Data\ExtracellularSummaries';
defaultSaveCond=true;
defaultSetMinSampThreshold=false;

%define optional inputs
addRequired(p,'DoseResponseCurve',@isstruct);
addRequired(p,'genotypes',@iscellstr);
addRequired(p,'doses',@isnumeric);
addOptional(p,'FigPath',defaultFigurePath,@ischar);
addOptional(p,'SaveCond',defaultSaveCond,@islogical);
addOptional(p,'MinSampThreshold',defaultSetMinSampThreshold,@islogical);

%parse
parse(p,DoseResponseCurve,genotypes,dose,varargin{:});

%reset defaults
figpath=p.Results.FigPath;
savecond=p.Results.SaveCond;
setminsamps=p.Results.MinSampThreshold;

%%

%initialize variable
removeidx=[];

%'dead' is defined as
deathspergenotype=nan(1,length(genotypes));
numpergenotype=nan(1,length(genotypes));
deadthresh=1/60; %ie one 'burst' per minute, this does not occur and would likely signify some minimal value of noise events
removeidx=[];
titlestr=('Samples per genotype: ');

for genotypeidx=1:length(genotypes)
    genotype=genotypes{genotypeidx};
    
    datavar=DoseResponseCurve.(genotype);
    
    header=datavar(1,:);
    tempindx=find(header=='doses');
    
    if ~isstring(dose)
        dose=string(dose);
    end
    doseindices=find(datavar(:,tempindx)==dose);
    
    tempindx=find(header=='avg_freq');
    tempvar=datavar(:,tempindx);
    allavgfreq=cellstr(tempvar(doseindices));
    allavgfreq_num=cellfun(@str2num,allavgfreq);
    
    deathbool=allavgfreq_num<deadthresh;
    
    numpergenotype(genotypeidx)=length(deathbool);
    
    numdeaths=sum(deathbool(:)==true);
    numlives=sum(deathbool(:)==false);
    
    deathrate=numdeaths/length(deathbool)*100;
    liverate=numlives/length(deathbool)*100;
    
    deathspergenotype(genotypeidx)=deathrate;
    livespergenotype(genotypeidx)=liverate;
    
    %define which indices to remove if there isn't at least 3 slices with this dose
    if length(deathbool)<3
        removeidx=[removeidx genotypeidx];
        titlestr=strcat(titlestr,sprintf(' %s=%d ',genotype,length(deathbool)));
    else
        titlestr=strcat(titlestr,sprintf(' %s=%d ',genotype,length(deathbool)));
    end
    
end

if ~setminsamps
    removeidx=[];
end

% remove these values from plotted inputs
genotypes(removeidx)=[];
livespergenotype(removeidx)=[];
deathspergenotype(removeidx)=[];
numpergenotype(removeidx)=[];

if ~sum(contains(genotypes,'wt'))
    return
elseif ~(length(genotypes)>2)
    return
end

h=figure;
plotlivespergenotype=flip(livespergenotype);
plotgenotypes=flip(genotypes);
barh(plotlivespergenotype,'FaceColor','w','EdgeColor','k','LineWidth',3);

%set graph params
set(gca,'Box','off')
set(gca,'FontName','Ariel');
set(gca,'FontSize',15);
set(gcf,'color','w');
set(gca,'yticklabel',plotgenotypes);
xlim([0 150])
xticks([0 50 100]);
xticklabs={'0%','50%','100%'};
set(gca,'xticklabel',xticklabs);
set(gca,'TickDir','out'); % The only other option is 'in'
title(titlestr)
ax=gca;
ax.TitleFontSizeMultiplier=0.5;
suptitlestr=sprintf('Rhythm Persists at %snM DAMGO',dose);
suptitle(suptitlestr)

if savecond
    cd(figpath);
    savefig(h,suptitlestr)
    saveas(h,suptitlestr,'pdf');
end

end
