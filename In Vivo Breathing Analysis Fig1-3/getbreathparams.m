function [MouseSummary,AllBreaths]= getbreathparams(animalID,varargin)
%this function segments breaths and gets relavent breathing
%parameters for all files in an animals folder, and organize them into
%accessible summary data structures.

%if pool object is open, close it so that it doesnt crash
poolobj=gcp('nocreate');
if ~isempty(poolobj)
    delete(poolobj)
end

%%create input parser

p=inputParser;

%define default conditions
defaultDataType='*.txt';
defaultSaveCond=true;
defaultSampleRate=1000;
defaultPoolSize=8;
defaultStartIndex=1;
defaultSampleLength=nan;
defaultAutoAdjust=false;
defaultDisplayTraces=false;

%parse optional inputs
addRequired(p,'animalID',@ischar);
addOptional(p,'SampleRate',defaultSampleRate,@isnumeric);
addOptional(p,'DataType',defaultDataType,@ischar);
addOptional(p,'SaveCond',defaultSaveCond,@islogical);
addOptional(p,'PoolSize',defaultPoolSize,@isnumeric);
addOptional(p,'StartIndex',defaultStartIndex,@isnumeric); %input defined in samples
addOptional(p,'SampleLength',defaultSampleLength,@isnumeric); %input defined in seconds
addOptional(p,'AutoCutArtifacts',defaultAutoAdjust,@islogical);
addOptional(p,'DisplayTraces',defaultDisplayTraces,@islogical);

%parse
parse(p,animalID,varargin{:});

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

%save parameters in summary file
MouseSummary.params=p.Results;

%% initiate parallel pool

parpool(poolsize);

%% get basic data info

%get data paths and conditions of experiment
[genotype,masterpathtodata,figpath,surgconds,expconds]=getmousepath(animalID);
cd(masterpathtodata);

% find all .txt (or other datatype) files in folder and subsequent subfolders
allFiles = dir(fullfile(masterpathtodata, '**',datatype));

disp(strcat('Initiating Analysis for: ',animalID));
fprintf('MouseSummary Save Condition is: %s\n',string(savecond))

%% iterates through each data file to segment breaths etc

tic
for i = 1 : length(allFiles)
    
    fName = allFiles(i).name;
    
    if ~strcmp(pwd,allFiles(i).folder) %if you're not in the right directory ie. its in a subfolder etc
        cd(allFiles(i).folder)
    end
    
    %get basic file information, this will depend on how you define your file names
    name=erase(fName,'.txt');
    mouseID = strsplit(name, '_');
    mouse = strcat(erase(mouseID{1},'cage'), '_', mouseID{2}); %gives cage and ID num
    surgcondition=string(mouseID{3});
    oxcondition=string(mouseID{4});
    drugcondition=string(mouseID{5});
    
    surgerycondition=string(expconds(lower(surgconds)==lower(surgcondition)));
    condition = strcat(surgerycondition, '_', oxcondition,'_',drugcondition); %before/after cre, normoxia v hypercapnia, morphine
    condition=strrep(condition,' ','_');
    toc
    disp(strcat("Analyzing... ",condition));
    tic
    
    %read in data
    A=readtable(fName,'ReadVariableNames',false);
    A=table2array(A);
    header=A(1:5,:);
    footer=A(end-2:end,:);
    A(1:5,:)=[];
    A(end-2:end,:)=[];
    sz=size(A);
    datalen=sz(1); % number of samples total
    
    %if the startindex is greater than the actual time sampled, let user
    %know and continue
    if startind>datalen
        disp('Data file does not fit expected dimentions...')
        disp('skipping file')
        %save data structures, in case this is last file
        savepath=strcat(masterpathtodata,'\',mouse,'_','summary','.mat');
        save(savepath,MouseSummary);
        savepath=strcat(masterpathtodata,'\',mouse,'_','raw','.mat');
        save(savepath,AllBreaths);
        continue
    end
    
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
    
    %% convert to usable data classes
    % do this by dividing in to chunks and running parallel conversion
    
    % split time/voltage into chunks
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
    
    % initialize new temporary variables
    newtime=NaN(numiter,25);
    newvoltage=NaN(numiter,25);
    
    % parallize conversion of datatypes, into new variables
    parfor iter=1:numiter
        newtime(iter,:)=converttoseconds(time(iter,:));
        newvoltage(iter,:)=cellfun(@str2double,voltage(iter,:));
    end
    
    % reshape and reset data vars
    time=reshape(newtime,[],1);
    voltage=reshape(newvoltage,[],1);
    voltage=voltage-mean(voltage); %center at 0 for ease
    clear newtime newvoltage
    
    % let me know if data does not appear as expected
    if time(1)~=.001 && startind>0
        disp('Data does not start at time 0:');
        disp(condition);
        time=time-time(1)+.001;
    end
    
    % let me know if data dies bit appear as expected
    if length(time)~= length(voltage)
        disp('Warning: error in data structure');
        %save data structures before break
        savepath=strcat(masterpathtodata,'\',mouse,'_','summary','.mat');
        save(savepath,MouseSummary);
        savepath=strcat(masterpathtodata,'\',mouse,'_','raw','.mat');
        save(savepath,AllBreaths);
        break
    end
    
    %% find major air flow artifacts
    
    airpress=movmean(voltage,2*fs); %take a moving average to see if airpressure is changing dramatically over time
    errorids=find(abs(airpress)>1.5)'; %should be centered at 0, this would be an extremely large event not related to breathing
    
    if ~isempty(errorids)
        %remove 1 sec around any error
        startremovallog=find(diff(errorids)>1);
        startremovallog=startremovallog+1;
        startremovals=[errorids(1) errorids(startremovallog)];
        stopremovallog=startremovallog-1;
        stopremovals=[errorids(stopremovallog) errorids(end)];
        numremovals=length(startremovals);
        removeintervals=[];
        for index=1:numremovals %remove data 1s before and after a mechanical error
            removeintervals=[removeintervals (startremovals(index)-1*fs):(stopremovals(index)+1*fs)]; %remove times +-1 second from event
        end
        removeintervals(find(removeintervals<1))=[]; %remove occassions that are impossible indices
        removeintervals(find(removeintervals>length(voltage)))=[]; %remove occassions that are impossible indices
        removeintervals=unique(removeintervals); %remove duplicates
        
        errorids=sort(removeintervals,'ascend');
        clear removeintervals
        
        %if you want, display these so you can determine if they're worth
        %dealing with
        if displaytrace
            rawdata=figure;
            
            plot(time,voltage,'k'); hold on;
            plot(time,airpress,'b');hold on;
            plot(time(errorids),ones(1,length(errorids)),'r.'); hold on;
            title(strrep(condition,'_',' '))
            xlabel('sec')
            ylabel('mV')
            
            if savecond
                savefig(rawdata,strcat(figpath,'\',condition,'_annotatedartifactremoval.fig'))
            end
            
        end
        close(rawdata)
        
        %if you want to allow cutting data based on air flow artifacts, adjust
        %accordingly
        if autoadjust
            %remove error instances from trace
            RemovedData.artifacts=voltage(errorids);
            voltage(errorids)=[];
            time(errorids)=[];
            time=eliminategaps(time,1/fs);
        end
    end
    
    clear airpress errorids numremovals
    
    %% identify initiation of breaths
    
    %define parameters for breath segmentation, these can be changed based on
    %recording system
    durThresh = 0.00002*10^5; %in seconds, minduration of inspiration
    inspAmpThresh = -0.5; %min inspiratory amplitude
    expAmpThresh = 0; %min expiratory amplitude
    
    %get breath start indices
    breathStarts = getbreathstarts(voltage, durThresh, inspAmpThresh, expAmpThresh);
    breathStarts=breathStarts(2:end); %sometimes first breath is cut off so not captured accurately, therefore throw it out
    
    %arrange all breaths in cell matrix, for ease of indexing
    breathmat=cell(1,length(breathStarts)-1);
    parfor breathindx=1:(length(breathStarts)-1)
        tempbreathstart=breathStarts(breathindx);
        tempbreathend=breathStarts(breathindx+1);
        tempbreath=voltage(tempbreathstart:tempbreathend)';
        breathmat{1,breathindx}=tempbreath;
    end
    
    %% use these to get respective breath parameters
    
    [inspPeak, expPeak, inspDur, expDur, inspVt, expVt, breathVt]=getbreathvals(breathmat);
    
    %% Remove events that might be classified as a 'breath' but are more likely some temporary noise fluctuation
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
    lowExpEvents=find(expPeak<1); %find events with small expiration
    
    %only include low amplitude events that also have significantly lower inspVt
    eliminateBreaths=find(abs(inspVt(lowAmpEvents))<abs(0.5*min(inspVt(highAmpEvents)))); %find which of these breaths also have a smaller VT than high amp events
    eliminateBreaths=lowAmpEvents(eliminateBreaths); %this is now indices of low amplitude, and low volume breaths
    
    %find which of these events also has small exp fluctuation
    reduction=ismember(eliminateBreaths,lowExpEvents); %this is indices of low amplitude, low volume, and low expiratory flow breaths, ie likely noise
    eliminateBreaths=eliminateBreaths(reduction);
    RemovedData.noiseflucts=breathStarts(eliminateBreaths); %save information for which breaths you are eliminating
    breathStarts(eliminateBreaths)=[];
    
    %if breaths were eliminated, recalculate parameters
    if ~isempty(eliminateBreaths)
        
        breathmat=cell(1,length(breathStarts)-1);
        parfor breathindx=1:(length(breathStarts)-1)
            tempbreathstart=breathStarts(breathindx);
            tempbreathend=breathStarts(breathindx+1);
            tempbreath=voltage(tempbreathstart:tempbreathend)';
            breathmat{1,breathindx}=tempbreath;
        end
        
        [inspPeak, expPeak, inspDur, expDur, inspVt, expVt, breathVt] = getbreathvals(breathmat);
        
    end
    
    %% validate segmentation by plotting waveforms and overlays of breathstart times over raw data
    
    if displaytrace %if you choose to
        
        exampledata=figure;
        
        %overlayed breaths
        subplot(1,2,1)
        breathindx=1;
        numbreaths=30;
        randbreathids=randsample(1:(length(breathStarts)-1),numbreaths);
        for breathindx=1:numbreaths
            subplot(1,2,1)
            plot(breathmat{1,randbreathids(breathindx)});
            hold on;
        end
        xlabel('ms');
        
        %annotated trace
        subplot(1,2,2)
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
        
        savefig(exampledata,strcat(figpath,'\',condition,'_breathexamples.fig'))
        
        close(exampledata)
        
    end
    
    
    %% save data into a usable data structure for each mouse
    %AllBreaths contains segmented breath matrices, seperated by recording
    %condition
    %MouseSummary contains all relevant breath parameters, seperated by
    %recording condition
    
    MouseSummary.(condition).inspPeak=inspPeak;
    MouseSummary.(condition).expPeak=expPeak;
    MouseSummary.(condition).inspDur=inspDur;
    MouseSummary.(condition).expDur=expDur;
    MouseSummary.(condition).Vt = inspVt;
    MouseSummary.(condition).expVt=expVt;
    MouseSummary.(condition).breathVt=breathVt;
    MouseSummary.(condition).minVt= inspVt*length(breathStarts)/(samplen/60); % Vt*numbreaths/min
    MouseSummary.(condition).Voltage = voltage; %I choose to save the entire voltage trace here as well, not really necessary and it does take up significant space but converting data structures also takes a long time
    MouseSummary.(condition).breathStartInd_indices=breathStarts; %in samples
    MouseSummary.(condition).breathStartInd =time(breathStarts);
    MouseSummary.(condition).datapath=name;
    MouseSummary.(condition).RemovedData=RemovedData;
    MouseSummary.(condition).emkainfo=header;
    AllBreaths.(condition)=breathmat;
    
    %save data structures
    if i==length(allFiles)
        savepath=strcat(masterpathtodata,'\',mouse,'_','summary','.mat');
        save(savepath,'MouseSummary');
        savepath=strcat(masterpathtodata,'\',mouse,'_','raw','.mat');
        save(savepath,'AllBreaths');
    end
    
    disp(sprintf('File %d/%d complete',i,length(allFiles))); %lets you know progress
    
end

end

function time= eliminategaps(time,expectedinterval)
for i=2:length(time)
    realinterval=time(i)-time(i-1);
    tol=1e-12;
    if ~(norm(realinterval-expectedinterval)<tol)
        time(i:end)=time(i:end)-(time(i)-time(i-1)-expectedinterval); %adjust so its continuous in time
    else
        continue
    end
end
end
