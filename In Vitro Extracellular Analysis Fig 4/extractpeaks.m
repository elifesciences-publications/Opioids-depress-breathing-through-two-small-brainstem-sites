function [ BurstSummary ] = extractpeaks(pathtodata,ch,fs,varargin)
% the purpose of this code is to take in extracellular data, and find peaks of activity

%% create input parser, parse inputs

p=inputParser;

% define default conditions
defaultManualCond=false; %do not need to manually annotate
defaultDataType='*.txt';
defaultSaveCond=true;
defaultFilterLowerBound=100;
defaultFilterUpperBound=4000;
defaultRawThreshold=nan;

% define optional inputs
addRequired(p,'pathtodata',@ischar);
addRequired(p,'ch',@ischar);
addRequired(p,'fs',@isnumeric);
addOptional(p,'AnalysisMode',defaultManualCond,@ischar);
addOptional(p,'DataType',defaultDataType,@ischar);
addOptional(p,'SaveCond',defaultSaveCond,@islogical);
addOptional(p,'FilterUpperBound',defaultFilterUpperBound,@isnumeric);
addOptional(p,'FilterLowerBound',defaultFilterLowerBound,@isnumber);
addOptional(p,'RawThreshold',defaultRawThreshold); %set a manual raw theshold

% parse
parse(p,pathtodata,ch,fs,varargin{:});

% reset defaults
if ~isempty(p.Results.AnalysisMode) && strcmp(p.Results.AnalysisMode,'manual')
    manualcond=true;
else
    manualcond=false;
end
datatype=p.Results.DataType;
savecond=p.Results.SaveCond;
absolute_Th_raw=p.Results.RawThreshold;

% initialize filter settings
upperbound=p.Results.FilterUpperBound;
lowerbound=p.Results.FilterLowerBound;
nyquist=fs/2;
[B,A]=butter(2,[lowerbound/nyquist upperbound/nyquist],'bandpass'); %2nd order bandpass butterworth filter

%% go to folder

cd(pathtodata)

%% read in all files in folder

[folder,experimentid]=fileparts(pathtodata);
datadir=dir(datatype);

%% alphabatize so all of your figs etc are in order

datadir=nestedSortStruct(datadir,{'name'},1); %ascending order ie reverse order
numconditions=length(datadir);

%% preload all files & preallocate burstSummary

for recid=1:numconditions
    if strcmp(datatype,'*.txt')
        
        txtfile=strcat(datadir(recid).folder,'/',datadir(recid).name);
        fileID = fopen(txtfile,'r');
        formatSpec = '%f';
        A=fscanf(fileID,formatSpec);
        
        % get rec condition name, this changes based on your naming scheme
        subexpid=datadir(recid).name;subexpid = erase(subexpid,".txt");
        temp=strsplit(subexpid,'_'); date=[temp{2} '/' temp{3} '/' temp{1}]; pclampid=temp{4};
        recCond=strjoin(temp(5:end),'_');
        
        DataVar.(recCond).raw=A;
        
    elseif strcmp(datatype,'*.mat')
        
        datafile=strcat(datadir(recid).folder,'/',datadir(recid).name);
        
        % get rec condition name, this changes based on your naming scheme
        subexpid=datadir(recid).name;subexpid = erase(subexpid,".mat");
        temp=strsplit(subexpid,'_'); date=[temp{2} '/' temp{3} '/' temp{1}]; pclampid=temp{4};
        recCond=strjoin(temp(5:end),'_');
        
        DataVar.(recCond).raw=load(datafile);
        
    end
end


%% convert these data into time and voltage signal

% for every recording session (20 min interval)
for recid=1:numconditions
    
    %get raw data & name of relevant file
    if strcmp(datatype,'*.txt')
        subexpid=datadir(recid).name;subexpid = erase(subexpid,".txt");
        temp=strsplit(subexpid,'_'); date=[temp{2} '/' temp{3} '/' temp{1}]; pclampid=temp{4};
        recCond=strjoin(temp(5:end),'_');
        
        A=DataVar.(recCond).raw;
        time=A(1:2:end)'; Vtrace=A(2:2:end)';
        Vtrace=Vtrace-mean(Vtrace); %normalize
        
        % remove first 10 seconds, bc pclamp generates artifacts at
        % beginning of recording
        time=time(10*fs:end);
        Vtrace=Vtrace(10*fs:end);
        
        %this may be case specific for our recordings
        if max(Vtrace)<0.05
            Vtrace=Vtrace*1000;
        end
        
        if contains(recCond,'baseline','IgnoreCase',true)
            rawmin=min(Vtrace)-0.1;
            rawmax=max(Vtrace)+0.1;
        end
        
        if isnan(absolute_Th_raw) & contains(recCond,'baseline','IgnoreCase',true) %if user has not manually set a threshold for raw data
            absolute_Th_raw=prctile(Vtrace,99.99); %set it to 99.9th percentile to be above SNR (bursts occur infrequently, and therefore are not large contributors to distribution of data)
        end
        
    elseif strcmp(datatype,'*.mat')
        subexpid=datadir(recid).name;subexpid = erase(subexpid,".mat");
        temp=strsplit(subexpid,'_'); date=[temp{2} '/' temp{3} '/' temp{1}]; pclampid=temp{4};
        recCond=strjoin(temp(5:end),'_');
        
        fullid=strcat(experimentid,{' '},subexpid);
        time=DataVar.(recCond).raw.c001_Time;
        Vtrace=DataVar.(recCond).raw.(ch);
        Vtrace=Vtrace-mean(Vtrace); %normalize
        
        % remove first 10 seconds, bc pclamp generates artifacts at
        % beginning of recording
        time=time(10*fs:end);
        Vtrace=Vtrace(10*fs:end);
        
        % this may be case specific for our recordings
        if max(Vtrace)<0.005
            Vtrace=Vtrace*1000;
        end
        
        if contains(recCond,'baseline','IgnoreCase',true)
            rawmin=min(Vtrace)-0.1;
            rawmax=max(Vtrace)+0.1;
        end
        
        if isnan(absolute_Th_raw) & contains(recCond,'baseline','IgnoreCase',true) %if user has not manually set a threshold for raw data
            absolute_Th_raw=prctile(Vtrace,99.99); %set it to 99.9th percentile to be above SNR (bursts occur infrequently, and therefore are not large contributors to distribution of data)
        end
        
    end
    
    %% process trace: bandpass, integral, baseline correction
    
    % filter data with simply bandpass filter, to get high freq data
    buttfilt=filtfilt(B,A,Vtrace);
    
    %integrate signal
    window=0.05*fs;
    newFs=fs/window;
    [newTime,intTrace]=movtrapz(window,time,abs(buttfilt));    %take moving integral at 50ms windows
    
    % remove filter width to eliminate filter artifacts
    removebuff=30*newFs;
    intTrace=intTrace(removebuff:end-removebuff);
    newTime=newTime(removebuff:end-removebuff);
    
    % baseline correct and center at 0
    intTrace=msbackadj(newTime',intTrace')';
    intTrace=intTrace-median(intTrace);
    
    if contains(recCond,'baseline','IgnoreCase',true)
        minint=min(intTrace)-0.5;
        maxint=max(intTrace)+0.5;
    end
    
    %% detect peaks on trace
    
    %optional: run peak finding on chunked data, if it is very variable
    %     chunklen=2*60*newFs;
    %     peakidx=runMassSpecWavelet_chunked(chunklen,intTrace,newTime,newFs,Fsmax,SNR_Th_filt);
    
    %run full data
    SNR_Th_filt=1.5;
    Fsmax=1.5; %lets say looking for roughly a signal at this freq
    [wCoefs,peakidx,peakInfo,removedPeaks]=MassSpecWavelet(intTrace,newTime,newFs,Fsmax,SNR_Th_filt,"mexh","automated");
    
    %define locs and pks
    locs=peakidx; % in samples
    pks=intTrace(locs);
    
    %% plot filtered data
    
    figure;
    f=gcf;
    xlabel('min');ylabel('mV');
    %     plot(time/60,Vtrace,'color',[0.5 0.5 0.5]);
    %     hold on;
    plot(newTime/60,intTrace,'k')
    hold on;
    
    %% go back to raw data, measure amplitude of bursts and burst time based
    % on max value +- 0.5 seconds, as well as do some mild thresholding on
    % raw trace
    
    %create zscore trace bc you also want threshold based on raw snr
    zVtrace=normalize_ib(Vtrace,1,'zscore');
    SNR_Th_raw=2;
    
    %find "real" burst time as max point of each burst
    window=0.5*fs;
    rawpeakindices=[]; %will be in samples
    rawpeakamplitudes=[];
    removeIdx=[];
    peakidx=1;
    while peakidx<=length(locs)
        t0=newTime(locs(peakidx)); %peak time from mass spec wavelet
        [minVal,peaktime]=min(abs(time-t0)); %find nearest value to peak time in raw trace
        index0=peaktime-window; %look at +- window around this time
        if index0<1 index0=1; end
        index1=peaktime+window;
        if index1>length(time) index1=length(time); end
        rawint=abs(Vtrace(index0:index1));
        
        maxpeak=find(rawint==max(rawint));%gives you position relative to all of payntfilt
        maxpeak=maxpeak(1);
        
        % optional removal of peaks based on raw data threshold (SNR and
        % absolute value)
        if zVtrace(index0+maxpeak)<SNR_Th_raw || rawint(maxpeak)<absolute_Th_raw
            removeIdx=[removeIdx peakidx];
        end
        
        %         if rawint(maxpeak)<absolute_Th_raw
        %             removeIdx=[removeIdx peakidx];
        %         end
        
        rawpeakindices=[rawpeakindices (index0+maxpeak)];
        rawpeakamplitudes=[rawpeakamplitudes rawint(maxpeak)];
        
        peakidx=peakidx+1;
        
    end
    
    %get actual peakamps and peaktimes
    rawpeakamptimes=time(rawpeakindices);
    
    %% remove extreme low threshold events in raw from raw, and filtered as well
    
    removeIdx=unique(removeIdx);
    removelen=length(removeIdx);
    if removelen>0
        tempstring=sprintf('Removed %d Low Threshold Events',removelen);
        disp(tempstring)
    end
    removedPeaks.orphanRawThresh=locs(removeIdx); %record what you have removed
    rawpeakindices(removeIdx)=[];
    rawpeakamplitudes(removeIdx)=[];
    rawpeakamptimes(removeIdx)=[];
    locs(removeIdx)=[];
    pks(removeIdx)=[];
    
    %redefine locs and pks as time, amp to save
    filtpeaktimes=newTime(locs);
    filtpeakamplitudes=pks;
    
    %% plot peaks on data to visualize
    
    set(0, 'currentfigure', f);  %# for figures
    plot(filtpeaktimes/60,filtpeakamplitudes,'k^','markerfacecolor',[0 0 1]);
    ylim([minint maxint]);
    hold on;
    
    
    %% save summary values in struct
    
    BurstSummary.(recCond).date=date;
    BurstSummary.(recCond).pclampid=pclampid;
    BurstSummary.(recCond).filtTrace=intTrace;
    BurstSummary.(recCond).filtTime=newTime;
    BurstSummary.(recCond).filtBurstTimes=filtpeaktimes;
    BurstSummary.(recCond).filtBurstAmplitudes=filtpeakamplitudes;
    BurstSummary.(recCond).rawBurstTimes=rawpeakamptimes;
    BurstSummary.(recCond).rawBurstAmplitudes=rawpeakamplitudes;
    BurstSummary.(recCond).removedPeaks=removedPeaks;
    
    %% save figures with peak overlays
    
    if savecond
        
        cd(folder)
        if ~exist(strcat(folder,'/','Figures'),'dir')
            mkdir('Figures')
        end
        cd('Figures')
        
        %go to channel folder
        figpath=pwd;
        if ~exist(strcat(figpath,'/',ch),'dir')
            mkdir(ch)
        end
        cd(ch)
        
        set(0, 'currentfigure', f);  %# for figures
        title(strcat(strrep(subexpid,'_',' '),': Annotated Peaks'));
        savefig(f,strcat('annotatedpeaks_',subexpid));
        saveas(f,strcat('annotatedpeaks_',subexpid,'.pdf'));
        
    end
    
    close all
    
    cd(pathtodata) %go back to main data path
    
end

%% save summary data

if savecond
    
    %go to mat files folder
    cd(folder)
    if ~exist(strcat(folder,'/','matFiles'),'dir')
        mkdir('matFiles')
    end
    cd('matFiles')
    
    %go to channel folder
    figpath=pwd;
    if ~exist(strcat(figpath,'/',ch),'dir')
        mkdir(ch)
    end
    cd(ch)
    
    save(strcat(folder,'/','matFiles','/',ch,'/','BurstSummary'),'BurstSummary')
    
    %write to text file the code conditions
    txtfilename=strcat(folder,'/','matFiles','/',ch,'/','BurstSummary','_params.txt');
    fileID=fopen(txtfilename,'w');
    fprintf(fileID,'Parameters for ExtracellullarXII_extractpeaks'); fprintf(fileID,'\n');
    fprintf(fileID,'Manual Annotation: %s',logicaltostring(manualcond));fprintf(fileID,'\n');
    fprintf(fileID,'Data Type Used: %s',datatype);fprintf(fileID,'\n');fprintf(fileID,'\n');
    fprintf(fileID,'Butterworth Upper Bound Frequency: %d',upperbound);fprintf(fileID,'\n');
    fprintf(fileID,'Filtered SNR Data Threshold: %d',SNR_Th_filt);fprintf(fileID,'\n');
    fprintf(fileID,'Raw SNR Data Threshold: %d',SNR_Th_raw);fprintf(fileID,'\n');
    fprintf(fileID,'Raw Vm Data Threshold: %d',absolute_Th_raw);fprintf(fileID,'\n');
    
end

end

function str=logicaltostring(logcond)
if logcond==(false)
    str='false';
elseif logcond==(true)
    str='true';
else
    str='error';
end
end

function newMatrix=normalize_ib(matrix, dim,varargin)

method=varargin(1);

matrix=abs(matrix);

if dim==1
elseif dim==2
    matrix=matrix';
end

sz=size(matrix); nrows=sz(1); ncols=sz(2);

newMatrix=zeros(nrows,ncols);
prctiles=[0.025 0.25 0.50 0.75 0.975]*100;
prctilevals=[0 prctiles];
for i=1:nrows
    stdev=std(matrix(i,:));
    meanval=mean(matrix(i,:));
    if strcmp(method,'normzscore')
        newMatrix(i,:)=(matrix(i,:)-meanval)./stdev;
        newMatrix(i,:)=matrix(i,:)./max(matrix(i,:));
    elseif strcmp(method,'prctile')
        prctile_array=prctile(matrix(i,:),prctiles);
        for j=2:length(prctile_array)
            newMatrix(i,find(matrix(i,:)>prctile_array(j-1)& matrix(i,:)<=prctile_array(j)))=prctilevals(j);
        end
        newMatrix(i,find(matrix(i,:)>prctile_array(end)))=prctilevals(end);
    elseif strcmp(method,'zscore')
        newMatrix(i,:)=(matrix(i,:)-meanval)./stdev; %generate z-score vals
    else
        error('this option is not available');
    end
end

if dim==1
elseif dim==2
    newMatrix=newMatrix';
end


end

function allpeakindices=runMassSpecWavelet_chunked(chunk,data,time,Fs,Fsmax,SNR_Th_filt)
datalen=length(data);
numchunks=floor(datalen/chunk);
overflow=int32(chunk*(datalen/chunk-floor(datalen/chunk)));
if overflow>0
    chunked_time=time(1:(end-overflow));
    chunked_data=data(1:(end-overflow));
    overflowTime=time((end-overflow+1):end);
    overflowData=data((end-overflow+1):end);
end
chunked_time=reshape(chunked_time,chunk,[])';
chunked_data=reshape(chunked_data,chunk,[])';

%GET PEAKS ALGORITHM
%do mexican hat convolution and either manual, or automated ridge
%finding
%per chunk.
SNR_Th_filt=1.0;
Fsmax=1.5; %lets say looking for roughly a signal at this freq
chunknum=1;
allpeakindices=[];
while chunknum<=numchunks
    subtime=chunked_time(chunknum,:);
    subtrace=chunked_data(chunknum,:);
    [wCoefs,peakidx,peakInfo,removedPeaks]=MassSpecWavelet(subtrace,subtime,Fs,Fsmax,SNR_Th_filt,"mexh","automated");
    
    peakidx=peakidx+((chunknum-1)*chunk);
    allpeakindices=[allpeakindices peakidx];
    chunknum=chunknum+1;
end
if overflow>0
    subtime=overflowTime;
    subtrace=overflowData;
    [wCoefs,peakidx,peakInfo,removedPeaks]=MassSpecWavelet(subtrace,subtime,Fs,Fsmax,SNR_Th_filt,"mexh","automated");
    
    peakidx=peakidx+(numchunks*chunk);
    allpeakindices=[allpeakindices peakidx];
end
end




