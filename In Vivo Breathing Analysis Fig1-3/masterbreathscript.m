% this master script executes a variety of analysis from 'Opioids depress
% breathing through two small brainstem sites' (Bachmutsky et al. 2020)

% specifically, this conducts analysis on in vivo breathing data recorded
%in a whole-body plethysmograph chamber, as seen in Figs. 1-3

% this code conducts analysis across mice listed in mouseids, can comment
%out sections you do not want to run, or run sections individually

% see doi: https://doi.org/10.1101/807297 for methods and relevant
% experimental parameters

%% plot single files if you need to validate any questionable recordings

txtfile='C:\Users\Iris Bachmutsky\Desktop\Bachmutskyetal.2020\In Vivo Example Breathing Data\155_2\Baseline\Hypercapnia\cage155_2_beforeinject_hypercapnia_saline.txt';
[BreathSummary,RawBreaths]=breathsegmentation_singlefile(txtfile,'SaveCond',true,'AutoCutArtifacts',true,'DisplayTraces',false);
plottrace_singlefile(txtfile,'SaveCond',true,'SavePath','C:\Users\Iris Bachmutsky\Desktop\Bachmutskyetal.2020\In Vivo Example Breathing Data');


%% choose which mice you want to run all analysis on

mouseids=mousegeneticsdirectory('all'); %ex: get all pbc ko, followed by kf ko mice

%% get the breath parameters for each mouse

for mouseindex=1:length(mouseids)
    mouseid=mouseids{mouseindex};
    [genotype,pathtodata,figpath,surgconds,expconds]=getmousepath(mouseid);
    cd(pathtodata);
    
    temp=strsplit(mouseid,'_');
    experimentorid=temp{1};
    mouseid2=strcat(temp{2},'_',temp{3});
    
    [~,AllBreaths]=getbreathparams(mouseid,'SaveCond',true,'AutoCutArtifacts',true,'DisplayTraces',true);
    
    %get these special expiratory params (pause)
    cd(pathtodata);
    splitexperations(mouseid2,AllBreaths,true)
    
end

%% plot breath histograms, ie histograms of different breath parameters,overlayed by drug condition
% see ex: Fig 1I

plotconditions={'exppause','instantaneous frequency','inspiratory peak','tidal volume'};
for mouseindex=1:length(mouseids)
    for plotindex=1:length(plotconditions);
        plotcondition=plotconditions{plotindex};
        mouseid=mouseids{mouseindex};
        
        plotbreathhistograms(mouseid,'PlotCondition',plotcondition)
        
    end
end

%% plot scatters of different parameters for each condition, assumes you want both saline and morphine overlayed
%see ex: Fig 1B/E

for mouseindex=1:length(mouseids)
    mouseid=mouseids{mouseindex};
    
    plotbreathscatters(mouseid,'XAxis','Instantaneous Frequency','YAxis','Peak Flow');
    plotbreathscatters(mouseid,'XAxis','Inspiratory Duration','YAxis','Peak Flow');
    plotbreathscatters(mouseid,'XAxis','Tidal Volume','YAxis','Peak Flow');
    
end

%% plot scatters of different parameters, but more customizable by condition and color
%see ex: Fig 2F

for mouseindex=1:length(mouseids)
    mouseid=mouseids{mouseindex};
    
    %order of plotconditions will define the order the images are stacked in
    
    %plot single ko overlay
    plotconditions={'Baseline_hypercapnia_saline','PBC_mOR_KO_hypercapnia_saline','Baseline_hypercapnia_morphine','PBC_mOR_KO_hypercapnia_morphine'};
    ColorMatrix={[0 0 0],[85,85,85]./255;[0 0 0],[85,85,85]./255;[150,0,0]./255,[255,0,0]./255;[0, 0, 255]./255,[77, 77, 255]./255};
    plotbreathscatters_preset(mouseid,'XAxis','Instantaneous Frequency','YAxis','Peak Flow','PlotConditions',plotconditions,'ColorMatrix',ColorMatrix,'SubFilename','PBC mOR KO');
    
    %plot double ko overlay
    plotconditions={'Baseline_hypercapnia_saline','PBC_mOR_KO_hypercapnia_saline','PBC_KF_mOR_KO_hypercapnia_saline','Baseline_hypercapnia_morphine','PBC_mOR_KO_hypercapnia_morphine','PBC_KF_mOR_KO_hypercapnia_morphine'};
    ColorMatrix={[0 0 0],[85,85,85]./255;[0 0 0],[85,85,85]./255;[0 0 0],[85,85,85]./255;[150,0,0]./255,[255,0,0]./255;[0, 0, 255]./255,[77, 77, 255]./255;[0, 179, 0]./255,[0, 255, 0]./255 };
    plotbreathscatters_preset(mouseid,'XAxis','Instantaneous Frequency','YAxis','Peak Flow','PlotConditions',plotconditions,'ColorMatrix',ColorMatrix,'SubFilename','PBC KF mOR KO');
    
    
end

%% plot percent time by breathing period
%see ex: Fig 1J

for mouseindex=1:length(mouseids)
    mouseid=mouseids{mouseindex};
    
    plotpercbreathingperiod(mouseid,'SaveCond',true)
    
end

%% plot overlays of breaths for binned conditions, to compare morphology

%pick condition you want to sort by, and which intervals to sort into
sortconditions={'inspiratory duration'};
inspdurs={[50 100],[100 150],[150 200]};

for mouseindex=1:length(mouseids)
    mouseid=mouseids{mouseindex};
    
    for sortindex=1:length(sortconditions);
        sortcondition=sortconditions{sortindex};
        
        switch sortcondition
            case 'inspiratory duration'
                bounds=inspdurs;
            case 'tidal volume'
                bounds=tidvolvals;
            case 'pause duration'
                bounds=pausedurvals;
        end
        
        for boundsindex=1:length(bounds)
            boundpair=bounds{boundsindex};
            lowbound=boundpair(1);
            upbound=boundpair(2);
            
            plotavgnoverlaybreaths(mouseid,'SaveCond',true,'SortCondition',sortcondition,'LowerBound',lowbound,'UpperBound',upbound,'PlotType','flow')
            
        end
    end
end

%% plot summary data, summarized by mouse across conditions
%see ex: Fig C,F,M

%set which conditions you want to plot, each will be saved in a seperate figure
plotconditions={'instfreq','peakflow','minVT','VT','exppause'};

% set plot conditions, ie exp. conditions you want to plot
pbckfkoexpconds={'Baseline','PBC mOR KO','PBC KF mOR KO'};

plotsummaryparams(mouseids,'PlotConditions',plotconditions,'SaveCond',true,'FolderName','PBC KF mOR KO','SurgConditions',pbckfkoexpconds);
