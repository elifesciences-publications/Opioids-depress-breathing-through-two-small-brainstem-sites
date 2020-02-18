function plottrace_singlefile(txtfile,varargin)
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
[filepath,name,ext]=fileparts(txtfile);
defaultSavePath=filepath;

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
addOptional(p,'SavePath',defaultSavePath,@ischar);

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
savepath=p.Results.SavePath;


%%
%get basic information, this is based on a personal data structure of
%choice from Bachmutsky et al. 2020

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
    time=time-time(1)+.001; % if you want to align with other measurment, this is not accurate
end

if length(time)~= length(voltage)
    disp('Warning: error in data structure');
    pause
end


%% plot breath data

    exampledata=figure;
  
    plot(time,voltage,'k');
    xlabel('sec');
    xlim([5*60 5.5*60]);
    ylim([-15 15]);
    
    suptitle(strrep(strcat(mouse,'_',condition),'_',' '))
    savestr=strcat(mouse,'_',condition,'_','breathtrace');
    
    if savecond 
        
    savefig(exampledata,strcat(savepath,'/',savestr,'.fig'));
    saveas(exampledata,strcat(savepath,'/',savestr,'.pdf'));

    end
    
end
