function [NumSlices]=plotdoseresponsecurve(DoseResponseCurve,genotypes,doses,varargin)
% the purpose of this code is to plot dose response curves of freq and
% amplitude across DAMGO doses

%%
%create input parser
p=inputParser;

%define default conditions
defaultPlotType='areacurve';
defaultFigurePath='C:\Users\Iris Bachmutsky\Desktop\Bachmutskyetal.2020\In Vitro Example Extracellular Data\ExtracellularSummaries';
defaultSaveCond=true;

%define optional inputs
addRequired(p,'DoseResponseCurve',@isstruct);
addRequired(p,'genotypes',@iscellstr);
addRequired(p,'doses',@isnumeric);
addOptional(p,'PlotType',defaultPlotType,@ischar);
addOptional(p,'FigPath',defaultFigurePath,@ischar);
addOptional(p,'SaveCond',defaultSaveCond,@islogical);

%parse
parse(p,DoseResponseCurve,genotypes,doses,varargin{:});

%reset defaults
plottype=p.Results.PlotType;
figpath=p.Results.FigPath;
savecond=p.Results.SaveCond;

%% initialize figure,vars

dosefig=figure;
[ha, pos] = tight_subplot(1,2,[.15 .15],[.15 .05],[.15 .2]);
numdoses=length(doses);
numgenotypes=length(genotypes);
maxfreq=[];
maxamp=[];
savestr=sprintf('DoseResponseCurve_%s',plottype);
tracecolors={'#0000FF';'#FF1234';'#565656';'#fea49f';'#fbaf08'};
tracecolors=hex2rgb(tracecolors);
plotvar=nan(1,numgenotypes);
stderror=@(x) std(x)/sqrt(length(x));
alphaarray=[0.2 0.5 0.7 0.9];
NumSlices=struct();


for genotypeidx=1:numgenotypes
    genotype=genotypes{genotypeidx};
    
    %%set up for plotting all avg freq, amplitudes
    
    %generate string for figure
    savestr=strcat(savestr,sprintf('_%s_filled',genotype));
    
    %get data for this genotype
    datavar=DoseResponseCurve.(genotype);
    header=datavar(1,:);
    
    %get all doses
    tempindx=find(header=='doses');
    dosevar=datavar(2:end,tempindx);
    dosevar=convertStringsToChars(dosevar);
    dosevar=cellfun(@str2num,dosevar);
    
    %get all avgfreq
    tempindx=find(header=='avg_freq');
    freqvar=datavar(2:end,tempindx);
    freqvar=convertStringsToChars(freqvar);
    freqvar=cellfun(@str2num,freqvar);
    
    %get all avgamp
    tempindx=find(header=='avg_amp');
    ampvar=datavar(2:end,tempindx);
    ampvar=convertStringsToChars(ampvar);
    ampvar=cellfun(@str2num,ampvar);
    
    %% arrange by freq/amp per dose that there are samples for
    
    %initialize vars
    
    if strcmp(plottype,'scatter')
        alldoses=[];
        allfreqs=[];
        allamps=[];
    else
        alldoses=nan(1,numdoses);
        allfreqs=cell(1,numdoses);
        allamps=cell(1,numdoses);
    end
    
    for doseidx=1:length(doses)
        dose=doses(doseidx);
        indices=find(dosevar==dose);
        
        if isempty(indices)
            continue;
        end
        
        if strcmp(plottype,'scatter')
            alldoses=[alldoses dosevar(indices)'];
            allfreqs=[allfreqs freqvar(indices)'];
            allamps=[allamps ampvar(indices)'];
        elseif strcmp(plottype,'areacurve')
            alldoses(doseidx)=dose;
            allfreqs{doseidx}=freqvar(indices);
            allamps{doseidx}=ampvar(indices);
        end
    end
    
    %% plot
    
    if strcmp(plottype,'scatter')
        %do a scatterplot version of dose response curve for freq and amp
        axes(ha(1))
        sigm_fit(alldoses,allfreqs,[],[],1); hold on;
        plotvar(genotypeidx)=scatter(alldoses,allfreqs,[],tracecolors(genotypeidx,:),'o','filled'); hold on;
        
        %reset axes
        if isempty(maxfreq)
            maxfreq=max(allfreqs);
            maxdose=max(alldoses);
        else
            tempvar=[allfreqs maxfreq];
            maxfreq=max(tempvar);
            tempvar=[alldoses maxdose];
            maxdose=max(tempvar);
        end
        ylim([0 maxfreq*1.2]);
        xlim([0 maxdose]);
        
        
        axes(ha(2))
        sigm_fit(alldoses,allamps,[],[],1); hold on;
        scatter(alldoses,allamps,[],tracecolors(genotypeidx,:),'o','filled'); hold on;
        %reset axes
        if isempty(maxamp)
            maxamp=max(allamps);
        else
            tempvar=[allamps maxamp];
            maxamp=max(tempvar);;
        end
        ylim([0 maxamp*1.2]);
        xlim([0 maxdose]);
        
    elseif strcmp(plottype,'areacurve')
        %do an area under the curve for average values with stdev
        
        num_slices=cellfun(@length,allfreqs);
        NumSlices.(genotype)=num_slices;
        
        freq_means=cellfun(@mean,allfreqs);
        freq_stderror=cellfun(@(x) stderror(x),allfreqs);
        amp_means=cellfun(@mean,allamps);
        amp_stderror=cellfun(@(x) stderror(x),allamps);
        
        axes(ha(1))
        plot(alldoses,freq_means,'color',tracecolors(genotypeidx,:)); hold on;
        e=errorbar(alldoses,freq_means,freq_stderror,freq_stderror,'o'); hold on;
        e.LineWidth=1;e.MarkerSize=5;e.Color=tracecolors(genotypeidx,:);e.MarkerEdgeColor=tracecolors(genotypeidx,:);e.MarkerFaceColor=tracecolors(genotypeidx,:);
        %fill in area under curve if wanted
        a=area(alldoses,freq_means,0); hold on;
        plotvar(genotypeidx)=a;
        a.FaceColor=tracecolors(genotypeidx,:);
        alphaval=alphaarray(genotypeidx);
        alpha(a,alphaval);
        xlim([0 500]);
        ylim([0 1.5]);
        
        axes(ha(2))
        plot(alldoses,amp_means,'color',tracecolors(genotypeidx,:)); hold on;
        f=errorbar(alldoses,amp_means,amp_stderror,amp_stderror,'o'); hold on;
        f.LineWidth=1.5;f.MarkerSize=5;f.Color=tracecolors(genotypeidx,:);f.MarkerEdgeColor=tracecolors(genotypeidx,:);f.MarkerFaceColor=tracecolors(genotypeidx,:);
        %fill in area under curve if wanted
        a=area(alldoses,amp_means,0);
        a.FaceColor=tracecolors(genotypeidx,:);
        alpha(a,alphaval);
        xlim([0 500]);
        ylim([0 1.5]);
        
    end
end


set(dosefig,'color','w');
axes(ha(1))

xlabel('nM DAMGO')
ylabel('Normalized Freq')
xticks(doses);
set(gca,'Box','off')
set(gca,'FontName','Ariel');
set(gca,'FontSize',10);
set(gca,'TickDir','out'); % The only other option is 'in'

% Construct a Legend with the data from the sub-plots
hL = legend(plotvar,genotypes);
newPosition = [0.825 0.8 0.1 0.05];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);

axes(ha(2))
xlabel('nM DAMGO')
ylabel('Normalized Amplitude')
xticks(doses);
set(gca,'Box','off')
set(gca,'FontName','Ariel');
set(gca,'FontSize',10);
set(gca,'TickDir','out'); % The only other option is 'in'

suptitle('Dose Response Curves')

if savecond
    cd(figpath);
    savefig(dosefig,savestr);
    saveas(dosefig,savestr,'pdf');
end

end
