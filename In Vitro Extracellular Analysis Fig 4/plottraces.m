function plottraces(pathtodata,ch,fs,varargin)
% the purpose of this code is to take in extracellular data, and plot
% integrated traces of voltage over time

%% create input parser, parse inputs

p=inputParser;

%define default conditions
defaultManualCond=false; %do not need to manually annotate
defaultDataType='*.txt';
defaultSaveCond=true;
defaultFilterLowerBound=100;
defaultFilterUpperBound=4000;
defaultPoolSize=8;

%define optional inputs
addRequired(p,'pathtodata',@ischar);
addRequired(p,'ch',@ischar);
addRequired(p,'fs',@isnumeric);
addOptional(p,'AnalysisMode',defaultManualCond,@ischar);
addOptional(p,'DataType',defaultDataType,@ischar);
addOptional(p,'SaveCond',defaultSaveCond,@islogical);
addOptional(p,'FilterUpperBound',defaultFilterUpperBound,@isnumeric);
addOptional(p,'FilterLowerBound',defaultFilterLowerBound,@isnumber);
addOptional(p,'PoolSize',defaultPoolSize,@isnumeric);

%parse
parse(p,pathtodata,ch,fs,varargin{:});

%reset defaults
if ~isempty(p.Results.AnalysisMode) && strcmp(p.Results.AnalysisMode,'manual')
    manualcond=true;
else
    manualcond=false;
end

datatype=p.Results.DataType;
savecond=p.Results.SaveCond;
poolsize=p.Results.PoolSize;

%get filter parameters
upperbound=p.Results.FilterUpperBound;
lowerbound=p.Results.FilterLowerBound;
nyquist=fs/2;
[B,A]=butter(2,[lowerbound/nyquist upperbound/nyquist],'bandpass'); %2nd order bandpass butterworth filter

%% go to folder

cd(pathtodata)

%% read in all files in folder

[folder,experimentid]=fileparts(pathtodata);
datadir=dir(datatype);

%%  alphabatize so all of your figs etc are in order

datadir=nestedSortStruct(datadir,{'name'},1); %ascending order ie reverse order
numconditions=length(datadir);

%% preload all files & preallocate burstSummary

for recid=1:numconditions
    
    if strcmp(datatype,'*.txt')
        
        txtfile=strcat(datadir(recid).folder,'/',datadir(recid).name);
        fileID = fopen(txtfile,'r');
        formatSpec = '%f';
        A=fscanf(fileID,formatSpec);
        
        %get rec condition name, this changes based on your naming scheme
        subexpid=datadir(recid).name;subexpid = erase(subexpid,".txt");
        temp=strsplit(subexpid,'_'); date=[temp{2} '/' temp{3} '/' temp{1}]; pclampid=temp{4};
        recCond=strjoin(temp(5:end),'_');
        
        DataVar.(recCond).raw=A;
        
    elseif strcmp(datatype,'*.mat')
        
        datafile=strcat(datadir(recid).folder,'/',datadir(recid).name);
        
        %get rec condition name, this changes based on your naming scheme
        subexpid=datadir(recid).name;subexpid = erase(subexpid,".mat");
        temp=strsplit(subexpid,'_'); date=[temp{2} '/' temp{3} '/' temp{1}]; pclampid=temp{4};
        recCond=strjoin(temp(5:end),'_');
        
        DataVar.(recCond).raw=load(datafile);
        
    end
end

%% convert these data into time and voltage signal

for recid=1:numconditions
    
    %get raw data & name of relevant file
    if strcmp(datatype,'*.txt')
        
        subexpid=datadir(recid).name;subexpid = erase(subexpid,".txt");
        disp(strcat("Initiate Plotting for: ",subexpid));
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
        
        rawmin=min(Vtrace)-0.1;
        rawmax=max(Vtrace)+0.1;
        
    elseif strcmp(datatype,'*.mat')
        subexpid=datadir(recid).name;subexpid = erase(subexpid,".mat");
        disp(strcat("Initiate Plotting for: ",subexpid));
        temp=strsplit(subexpid,'_'); date=[temp{2} '/' temp{3} '/' temp{1}]; pclampid=temp{4};
        recCond=strjoin(temp(5:end),'_');
        
        time=DataVar.(recCond).raw.c001_Time;
        Vtrace=DataVar.(recCond).raw.(ch);
        Vtrace=Vtrace-mean(Vtrace); %normalize
        
        % remove first 10 seconds, bc pclamp generates artifacts at
        % beginning of recording
        time=time(10*fs:end);
        Vtrace=Vtrace(10*fs:end);
        
        %this may be case specific for our recordings
        if max(Vtrace)<0.005
            Vtrace=Vtrace*1000;
        end
        
        rawmin=min(Vtrace)-0.1;
        rawmax=max(Vtrace)+0.1;
        
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
    newTime=newTime-(newTime(1));
    
    %baseline correct and center at 0
    intTrace=msbackadj(newTime',intTrace')';
    intTrace=intTrace-median(intTrace);
    
    %% plot integrated,filtered data
    
    figure;
    g=gcf;
    plot(newTime,intTrace,'color','k');
    title(strcat(strrep(subexpid,'_',' '),': Integrated Trace'));
    xlabel('sec');ylabel('mV');
    ylim([intmin intmax]);
    hold on;
    
    %% plot raw data
    
    figure;
    f=gcf;
    plot(time,Vtrace,'color','k');
    title(strcat(strrep(subexpid,'_',' '),': Raw Trace'));
    xlabel('sec');ylabel('mV');
    ylim([rawmin rawmax]);
    hold on;
    
    %% save figures
    
    %go to figures folder
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
    
    savefig(g,strcat('integrated_',subexpid));
    saveas(g,strcat('integrated_',subexpid,'.pdf'));
    close(g)
    savefig(f,strcat('raw_',subexpid));
    saveas(f,strcat('raw_',subexpid,'.pdf'));
    close(f)
    
    cd(pathtodata) %go back to main data path
end


end


