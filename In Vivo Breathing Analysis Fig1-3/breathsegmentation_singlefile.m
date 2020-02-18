function [BreathSummary,RawBreaths]= breathsegmentation_singlefile(txtfile,varargin)
%conduct breath segmentation on single recording

%%
%create input parser
p=inputParser;
%define default conditions
defaultDataType='*.txt';
defaultSaveCond=true;
defaultSampleRate=1000;
defaultPoolSize=8;
defaultStartIndex=1; %for where you want to start analyzing relative to onset of recording. in samples.
defaultSampleLength=nan;
defaultAutoAdjust=false; %if you want to automatically adjust recording artifacts
defaultDisplayTraces=false; %if you want to display validation of breath segmentation

%parse optional inputs
addRequired(p,'animalID',@ischar);
addOptional(p,'SampleRate',defaultSampleRate,@isnumeric);
addOptional(p,'DataType',defaultDataType,@ischar);
addOptional(p,'SaveCond',defaultSaveCond,@islogical);
addOptional(p,'PoolSize',defaultPoolSize,@isnumeric);
addOptional(p,'StartIndex',defaultStartIndex,@isnumeric);
addOptional(p,'SampleLength',defaultSampleLength,@isnumeric); %input defined in seconds
addOptional(p,'AutoCutArtifacts',defaultAutoAdjust,@islogical);
addOptional(p,'DisplayTraces',defaultDisplayTraces,@islogical);

%parse
parse(p,txtfile,varargin{:});

%reset defaults
datatype=p.Results.DataType;
savecond=p.Results.SaveCond;
fs=p.Results.SampleRate;
poolsize=p.Results.PoolSize;
startind=p.Results.StartIndex;
subsamp=p.Results.SampleLength;
subsamp=subsamp*fs; % input now defined in num samples
autoadjust=p.Results.AutoCutArtifacts;
displaytrace=p.Results.DisplayTraces;


%%
%get basic information, this is based on a personal data structure of
%choice from Bachmutsky et al. 2020

[filepath,name,ext]=fileparts(txtfile);
mouseID = strsplit(name, '_');
mouse = strcat(mouseID{1}, '_', mouseID{2}); %gives cage and ID num
surgcondition=string(mouseID{3}); %which injection number it is
oxcondition=string(mouseID{4}); %whether it was recorded in hypercapnia or normoxia
drugcondition=string(mouseID{5}); %which drug was aministered, for ex: saline, or morphine
condition = strcat(surgcondition, '_', oxcondition,'_',drugcondition); %before/after cre, normoxia v hypercapnia, morphine
condition=strrep(condition,' ','_');

%%
%read in data, convert to appropriate data type

A=readtable(txtfile,'ReadVariableNames',false);
A=table2array(A);
header=A(1:5,:);
footer=A(end-2:end,:);
A(1:5,:)=[];
A(end-2:end,:)=[];
sz=size(A);
datalen=sz(1); % number of samples total

%extract time and voltage information from table
if isnan(subsamp)
    time=A(startind:end,1);
    voltage=A(startind:end,2);
    samplen=(datalen-startind+1)/fs; %this is in seconds
else
    try
        time=A(startind:(startind+subsamp),1);
        voltage=A(startind:(startind+subsamp),2);
        samplen=(subsamp-startind+1)/fs;
    catch
        disp('Data file does not fit expected dimentions...')
        disp(condition);
        %take all available data instead, mark that you did so
        time=A(startind:end,1);
        voltage=A(startind:end,2);
        samplen=(datalen-startind+1)/fs;
        disp(strcat('new sample length:', num2str(samplen)));
    end
end
clear A

%with parallelization, convert from char arrays to double
newdatalen=length(time);
overflow=int32(25*(newdatalen/25-floor(newdatalen/25)));
if overflow>0
    time=time(1:(end-overflow));
    voltage=voltage(1:(end-overflow));
end
time=reshape(time,[],25);
voltage=reshape(voltage,[],25);
sz=size(time);
numiter=sz(1);
newtime=NaN(numiter,25);
newvoltage=NaN(numiter,25);
parfor iter=1:numiter
    newtime(iter,:)=converttoseconds(time(iter,:));
    newvoltage(iter,:)=cellfun(@str2double,voltage(iter,:));
end
time=reshape(newtime,[],1);
voltage=reshape(newvoltage,[],1);
voltage=voltage-mean(voltage); %center at 0 for ease
clear newtime newvoltage

if time(1)~=.001 && startind>0
    disp('Data does not start at time 0:');
    disp(condition);
    time=time-time(1)+.001;
end

if length(time)~= length(voltage)
    disp('Warning: error in data structure');
    pause
end


%%
%identify initiation of breaths

%define parameters for breath segmentation, these can be changed based on
%recording system
durThresh = 0.00002*10^5; %in seconds, minduration of inspiration
inspAmpThresh = -0.5; %min inspiratory amplitude
expAmpThresh = 0; %min expiratory amplitude

%get breath start indices
breathStarts = getbreathstarts(voltage, durThresh, inspAmpThresh, expAmpThresh);
breathStarts=breathStarts(2:end); %sometimes first breath is cut off so not captured accurately, therefore throw it out

%arrange all breaths in cell matrix, for ease of indexing
breathmat=cell(1,length(breathStarts));
parfor breathindx=1:(length(breathStarts)-1)
    tempbreathstart=breathStarts(breathindx);
    tempbreathend=breathStarts(breathindx+1);
    tempbreath=voltage(tempbreathstart:tempbreathend)';
    breathmat{1,breathindx}=tempbreath;
end


%%
%use these to get respective breath parameters

[inspPeak, expPeak, inspDur, expDur, inspVt, expVt, breathVt]=getbreathvals(breathmat);

%%
%Remove events that might be classified as a 'breath' but are more likely some temporary noise fluctuation
%Specifically, remove breaths that are
% 1) low amplitude
% 2) low volume
%and
% 3) low/no expiration

%divide breaths into low and high amplitude
newInspAmpThresh=2*inspAmpThresh;
lowAmpEvents=find(inspPeak>newInspAmpThresh); %find events that are smaller (now not making a threshold of -1)
highAmpEvents=find(inspPeak<=newInspAmpThresh);

%find all breaths with low amplitude expiration
lowExpEvents=find(expPeak<1); %if there is a very small inspiration but a big expiration... i'll save it regardless

%only include low amplitude events that also have significantly lower inspVt
eliminateBreaths=find(abs(inspVt(lowAmpEvents))<abs(0.25*min(inspVt(highAmpEvents)))); %find which of these breaths also have a smaller VT than high amp events
eliminateBreaths=lowAmpEvents(eliminateBreaths);

%find which of these events also has small exp fluctuation
reduction=ismember(eliminateBreaths,lowExpEvents);
eliminateBreaths=eliminateBreaths(reduction);
RemovedData.noiseflucts=breathStarts(eliminateBreaths); %save information for which breaths you are eliminating
breathStarts(eliminateBreaths)=[];

%if breaths were eliminated, recalculate parameters
if ~isempty(eliminateBreaths)
    
    breathmat=cell(1,length(breathStarts));
    parfor breathindx=1:(length(breathStarts)-1)
        tempbreathstart=breathStarts(breathindx);
        tempbreathend=breathStarts(breathindx+1);
        tempbreath=voltage(tempbreathstart:tempbreathend)';
        breathmat{1,breathindx}=tempbreath;
    end
    
    [inspPeak, expPeak, inspDur, expDur, inspVt, expVt, breathVt] = getbreathvals(breathmat);
    
end

%%
%validate segmentation by plotting waveforms and overlays of
%breathstart times over raw data

if displaytrace %if you choose to
    
    exampledata=figure;
    
    %overlayed breaths
    subplot(1,2,1)
    breathindx=1;
    numbreaths=30;
    randbreathids=randsample(1:length(breathStarts),numbreaths);
    while breathindx<numbreaths
        subplot(1,2,1)
        plot(breathmat{1,randbreathids(breathindx)});
        hold on;
        breathindx=breathindx+1;
        if breathindx>=length(breathStarts)
            breathindx=100000; %end while loop
        end
    end
    xlabel('ms');
    
    %annotated trace
    subplot(1,2,2)
    %     plot(oldtime,oldvoltage,'b'); hold on;
    plot(time(1:60*fs),voltage(1:60*fs),'k'); hold on;
    numbreathstarts=sum(breathStarts<=60*fs);
    for index=1:numbreathstarts
        plot([time(breathStarts(index)),time(breathStarts(index))],[-50,+50],'--.r');
        hold on;
    end
    xlim([time(1),time(5*fs)]);
    ylim([-15 15]);
    xlabel('sec');
    suptitle(strrep(strcat(mouse,'_',condition),'_',' '))
    
end

%%
%output into usable data structure

    BreathSummary.inspPeak=inspPeak;
    BreathSummary.expPeak=expPeak;
    BreathSummary.inspDur=inspDur;
    BreathSummary.expDur=expDur;
    BreathSummary.Vt = inspVt;
    BreathSummary.expVt=expVt;
    BreathSummary.breathVt=breathVt;
    BreathSummary.minVt= inspVt*length(breathStarts)/(samplen/60); % Vt*numbreaths/min
    BreathSummary.Voltage = voltage; %I choose to save the entire voltage trace here as well, not really necessary and it does take up significant space but converting data structures also takes a long time
    BreathSummary.breathStartInd_indices=breathStarts; %in samples
    BreathSummary.breathStartInd =time(breathStarts);
    BreathSummary.datapath=txtfile;
    BreathSummary.RemovedData=RemovedData;
    BreathSummary.emkainfo=header;
    RawBreaths=breathmat;
    
end
