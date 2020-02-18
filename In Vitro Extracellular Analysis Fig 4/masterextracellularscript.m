
% this master script executes a variety of analysis from 'Opioids depress
% breathing through two small brainstem sites' (Bachmutsky et al. 2020)

% specifically, this conducts analysis on extracellular slice recordings of
% the PreBotzinger rhythm, as seen in Fig. 4

% this code conducts analysis across recordings listed in recids, can
% comment out sections you do not want to run, or run sections individually

% see doi: https://doi.org/10.1101/807297 for methods and relevant
% experimental parameters


%% plot single files if you need to visualize any recordings before conducting analysis

datafile='C:\Users\Iris Bachmutsky\Desktop\Bachmutskyetal.2020\In Vitro Example Extracellular Data\rec38_2019_04_16_vglut2_wt_oprflox_DAMGO\Raw Data\2019_04_16_0004_DAMGO_50nM.mat';
ch='c002_IN_0'; % name of channel recorded on pclamp
fs=10000; % sampling frequency
plotsingletrace(datafile,ch,fs);

%% choose which recordings you want to conduct analysis on

recids=recordinggeneticsdirectory('all'); %ex: get all vglut2-cre, opr fl/fl slice recordings

%% get all relevant peak information, plot traces.
% see ex: Fig 4E

numrecs=length(recids);
for rec=1:numrecs
    recid=recids{rec};
    [pathtodata,figpath,~,datatype,inputchs,qualchs,genotypes,filtered]=getextracellular(recid);
    fs=10000;

    numchans=length(qualchs);
    if numchans~=length(genotypes)
        error('data directory not properly updated')
    end

    for chan=1:numchans
        chid=qualchs{chan};

        plottraces(pathtodata,chid,fs,'DataType',datatype,'SaveCond',true);
        extractpeaks(pathtodata,chid,fs,'DataType',datatype,'SaveCond',true)
    end
end

%% plot poincare plots of instantaneous frequency

numrecs=length(recids);
for rec=1:numrecs
    recid=recids{rec};
    fs=10000;

    [~,~,~,datatype,inputchs,qualchs,genotypes,~]=getextracellular(recid);


    numchans=length(qualchs);
    if numchans~=length(genotypes)
        error('data directory not properly updated')
    end

    for chan=1:numchans
        chid=qualchs{chan};
        genotype=genotypes{chan};

        plotpoincare(recid,chid,fs,'DataType',datatype,'Genotype',genotype);
    end
end

%% get relative dose responses for all conditions across genotypes (using
% last five minutes of recording data, to ensure full saturation of drug)

summarydatapath='C:\Users\Iris Bachmutsky\Desktop\Bachmutskyetal.2020\In Vitro Example Extracellular Data\ExtracellularSummaries';
if ~exist(summarydatapath)
    mkdir(summarydatapath)
end

%get raw values
DoseResponseCurve=struct();
numrecs=length(recids);
for rec=1:numrecs
    recid=recids{rec};
    fs=10000;
    
    [~,~,~,datatype,inputchs,qualchs,genotypes,~]=getextracellular(recid);
    
    
    numchans=length(qualchs);
    if numchans~=length(genotypes)
        error('data directory not properly updated')
    end
    
    for chan=1:numchans
        chid=qualchs{chan};
        genotype=genotypes{chan};
        
        [doses,avg_freq,avg_instfreq,avg_amp]=getdoseresponsecurve(recid,chid,fs,'Normalize','none');
        
        header=["recid","chid","doses","avg_freq","avg_instfreq","avg_amp"];
        doselen=length(doses);
        tempvar=strings(doselen,6);
        tempvar(1:doselen,1)=repmat(recid,doselen,1);
        tempvar(1:doselen,2)=repmat(chid,doselen,1);
        tempvar(1:doselen,3)=doses';
        tempvar(1:doselen,4)=avg_freq';
        tempvar(1:doselen,5)=avg_instfreq';
        tempvar(1:doselen,6)=avg_amp';
        
        %add header to this struct if it doesn't exist
        if ~isfield(DoseResponseCurve,genotype)
            DoseResponseCurve.(genotype)=header;
        end
        
        DoseResponseCurve.(genotype)=[DoseResponseCurve.(genotype); tempvar];
        
    end
end
cd(summarydatapath);
save('DoseResponseCurve_raw.mat','DoseResponseCurve')

%get values normalized to baseline
DoseResponseCurve=struct();
numrecs=length(recids);
for rec=1:numrecs
    recid=recids{rec};
    fs=10000;
    
    [~,~,~,datatype,inputchs,qualchs,genotypes,~]=getextracellular(recid);
    
    
    numchans=length(qualchs);
    if numchans~=length(genotypes)
        error('data directory not properly updated')
    end
    
    for chan=1:numchans
        chid=qualchs{chan};
        genotype=genotypes{chan};
        
        [doses,avg_freq,avg_instfreq,avg_amp]=getdoseresponsecurve(recid,chid,fs,'Normalize','baseline');
        
        header=["recid","chid","doses","avg_freq","avg_instfreq","avg_amp"];
        doselen=length(doses);
        tempvar=strings(doselen,6);
        tempvar(1:doselen,1)=repmat(recid,doselen,1);
        tempvar(1:doselen,2)=repmat(chid,doselen,1);
        tempvar(1:doselen,3)=doses';
        tempvar(1:doselen,4)=avg_freq';
        tempvar(1:doselen,5)=avg_instfreq';
        tempvar(1:doselen,6)=avg_amp';
        
        %add header to this struct if it doesn't exist
        if ~isfield(DoseResponseCurve,genotype)
            DoseResponseCurve.(genotype)=header;
        end
        
        DoseResponseCurve.(genotype)=[DoseResponseCurve.(genotype); tempvar];
        
    end
end
cd(summarydatapath);
save('DoseResponseCurve_normtobaseline.mat','DoseResponseCurve')

%% plot % dead (or alive) for each genotype
% ie. how many slices have an active rhythm at a dose of DAMGO
% see ex: Fig 4F

summarydatapath='C:\Users\Iris Bachmutsky\Desktop\Bachmutskyetal.2020\In Vitro Example Extracellular Data\ExtracellularSummaries';
if ~exist(summarydatapath)
    mkdir(summarydatapath)
end
genotypes={'wt','gad2','vglut2','dbx1'};%allgenotypes;

cd(summarydatapath);
temp=load('DoseResponseCurve_raw.mat');
DoseResponseCurve=temp.DoseResponseCurve;
doses=[50]; %nM DAMGO
for doseidx=1:length(doses)
    dose=doses(doseidx);
    
    plotpercalive(DoseResponseCurve,genotypes,dose);
    
end

%% plot dose response curve for each genotype
% see ex: Fig 4J

summarydatapath='C:\Users\Iris Bachmutsky\Desktop\Bachmutskyetal.2020\In Vitro Example Extracellular Data\ExtracellularSummaries';
if ~exist(summarydatapath)
    mkdir(summarydatapath)
end
cd(summarydatapath);
load('DoseResponseCurve_normtobaseline.mat');

genotypes={'vglut2','dbx1','wt'};
doses=[0 50 100 500];
[numslices]=plotdoseresponsecurve(DoseResponseCurve,genotypes,doses); 
