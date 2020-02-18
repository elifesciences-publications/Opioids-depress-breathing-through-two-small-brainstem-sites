function extracellular_plotpoincare(recid,ch,fs,varargin)
% plot poincare plots of frequency for each dose

%%
% create input parser
p=inputParser;

%define default conditions
defaultDataType='*.txt';
defaultSaveCond=true;
defaultPoolSize=8;
defaultGenotype='wt';

%define optional inputs
addRequired(p,'recid',@ischar);
addRequired(p,'ch',@ischar);
addRequired(p,'fs',@isnumeric);
addOptional(p,'DataType',defaultDataType,@ischar);
addOptional(p,'Genotype',defaultGenotype,@ischar);
addOptional(p,'SaveCond',defaultSaveCond,@islogical);
addOptional(p,'PoolSize',defaultPoolSize,@isnumeric);

%parse
parse(p,recid,ch,fs,varargin{:});

%reset defaults
datatype=p.Results.DataType;
savecond=p.Results.SaveCond;
poolsize=p.Results.PoolSize;
genotype=p.Results.Genotype;


%%
%get summary data for recording
[pathtodata,figpath,sumdatapath,~,~,~,~]=getextracellular(recid);

cd(sumdatapath);

cd(ch);

load('BurstSummary.mat');


%% get data for each condition

conditions=fieldnames(BurstSummary);
condlen=length(conditions);

%generate figure
g=figure;

%check how many recs have enough data to plot
plotlen=0;
for condidx=1:condlen
    
    %get name of rec
    condition=string(conditions(condidx));
    
    %get data for that recording
    datastruct=BurstSummary.(condition);
    
    %get burst times
    bursttime=datastruct.rawBurstTimes;
    
    if length(bursttime)<2
        continue
    end
    
    plotlen=plotlen+1;
    
end

plotlen=plotlen+1;

%% plot data for each condition 

%determine plotting parameters
nrows=ceil(plotlen/4);
collen=4;

%initialize vars
binedges=0:0.01:1;
norminstanfreq_summary=[];
nplus1_summary=[];
plotidx=0;
[ha, pos] = tight_subplot(nrows,collen);

for condidx=1:condlen
    
    %get name of rec
    condition=string(conditions(condidx));
    
    %get data for that recording
    datastruct=BurstSummary.(condition);
    
    %get burst times
    bursttime=datastruct.rawBurstTimes;
    
    if length(bursttime)<2
        continue
    end
    
    %get instan freq (Hz)
    instanfreq=1./diff(bursttime);
    norminstanfreq=instanfreq/max(instanfreq);
    norminstanfreq_summary=[norminstanfreq_summary norminstanfreq];
    
    %generate y axis of poincare plot
    nplus1=[0 norminstanfreq(1:end-1)];
    nplus1_summary=[nplus1_summary nplus1];
    
    %generate histogram for plotting
    poinc=histcounts2(nplus1,norminstanfreq,binedges,binedges);
    
    %plot
    plotidx=plotidx+1;
    axes(ha(plotidx));
    imagesc(binedges,binedges,poinc);
    axis equal;
    set(gca,'Xlim',binedges([1 end]),'Ylim',binedges([1 end]),'YDir','normal');
    colormap jet
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    title(condition,'Fontsize',8)
    
end

%generate summary plot of poincare data
poinc=histcounts2(nplus1_summary,norminstanfreq_summary,binedges,binedges);

%plot
plotidx=plotidx+1;
axes(ha(plotidx));
imagesc(binedges,binedges,poinc);
axis equal;
set(gca,'Xlim',binedges([1 end]),'Ylim',binedges([1 end]),'YDir','normal');
colormap jet
set(gca,'xtick',[])
set(gca,'ytick',[])
title('Overlay','Fontsize',8)
suptitle("Poincare Plots for Instan Freq")

%save figure
if savecond
    cd(figpath)
    cd(ch)
    savefig(g,'Poincare Plots for Instan Freq');
end

close all

end