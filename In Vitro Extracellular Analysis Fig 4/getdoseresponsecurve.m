function [alldoses,allavgfreq,allavginstanfreq,allavgamp]=extracellularXII_getdoseresponsecurve(recid,ch,fs,varargin)
%the purpose of this function is to get the freq/amplitude values for each dose

%%
%create input parser
p=inputParser;

%define default conditions
defaultNorm='baseline';
defaultSaveCond=true;
defaultPoolSize=8;

%define optional inputs
addRequired(p,'recid',@ischar);
addRequired(p,'ch',@ischar);
addRequired(p,'fs',@isnumeric);
addOptional(p,'SaveCond',defaultSaveCond,@islogical);
addOptional(p,'PoolSize',defaultPoolSize,@isnumeric);
addOptional(p,'Normalize',defaultNorm,@ischar);

%parse
parse(p,recid,ch,fs,varargin{:})

%reset defaults
normalization=p.Results.Normalize;

%% initialize variables

[pathtodata,figpath,sumdatapath,~,~,~,~]=getextracellular(recid);
cd(sumdatapath);
cd(ch);
load('BurstSummary.mat');
fnames=fieldnames(BurstSummary);
alldoses=nan(1,length(fnames));
allavgfreq=nan(1,length(fnames));
allavgamp=nan(1,length(fnames));
allavginstanfreq=nan(1,length(fnames));
dur=60*5; %five minutes, duration im extracting data from

%% get values to normalize recording features to (ie. normalize to baseline freq/amp, recovery or none)

if strcmpi(normalization,'baseline')
    
    %get baseline freq,amp for normalization
    isField=contains(fnames,'baseline','IgnoreCase',true);
    isField=find(isField,true,'last');
    
    %get the data
    if any(isField)
        datastruct= BurstSummary.(fnames{isField});
    else
        datastruct = [];
    end
    
    %tell me if data is missing
    if isempty(datastruct) %if still not found display that it is not present
        disp('missing baseline') 
        return
    end
    
    %extract avg values for instan freq, freq, and amp for last dur min of
    %recording
    [avg_instan_freq,avg_freq,avg_instan_amp]=getrecdata(datastruct,dur);
    
    norm_avgfreq=avg_freq;
    norm_avginstanfreq=avg_instan_freq;
    norm_avgamp=avg_instan_amp;
    
elseif strcmpi(normalization,'recovery')
    %get naloxone recover freq, amp for normalization
    isField=contains(fnames,'naloxone','IgnoreCase',true);
    isField=find(isField,true,'last');
    
    %get the data
    if any(isField)
        datastruct= BurstSummary.(fnames{isField});
    else
        isField=contains(fnames,'recovery','IgnoreCase',true);
        isField=find(isField,true,'last');
        if any(isField)
            datastruct= BurstSummary.(fnames{isField});
        else
            datastruct = [];
        end
    end
    %tell me if data is missing
    if isempty(datastruct) %if still not found display that it is not present
        disp('missing naloxone recovery')
        return
    end
    
    [avg_instan_freq,avg_freq,avg_instan_amp]=getrecdata(datastruct,dur);
    
    norm_avgfreq=avg_freq;
    norm_avginstanfreq=avg_instan_freq;
    norm_avgamp=avg_instan_amp;
    
elseif strcmp(normalization,'none')
    
    norm_avgfreq=1;
    norm_avginstanfreq=1;
    norm_avgamp=1;
    
else
    
    error('incorrect input for normalization')
    
end

%% iterate over dosages/fields per rec and normalize data

for doseindex=1:length(fnames)
    doseid=fnames(doseindex);
    doseid=doseid{:};
    
    % no need to analyze recovery condition
    if contains(doseid,'naloxone','IgnoreCase',true) || contains(doseid,'recover','IgnoreCase',true)
        continue
    end
    
    %get nM DAMGO for recording
    dose=getdose(doseid);
    alldoses(doseindex)=dose;
    
    %get data for that recording
    datastruct=BurstSummary.(doseid);
    
    [avg_instan_freq,avg_freq,avg_instan_amp]=getrecdata(datastruct,dur);
    
    allavgfreq(doseindex)=avg_freq/norm_avgfreq;
    allavgamp(doseindex)=avg_instan_amp/norm_avgamp;
    allavginstanfreq(doseindex)=avg_instan_freq/norm_avginstanfreq;
    
end

%remove nans from fields bc sometimes there are multiple nalox/recovery
%files
alldoses=alldoses(~isnan(alldoses));
allavgfreq=allavgfreq(~isnan(allavgfreq));
allavginstanfreq=allavginstanfreq(~isnan(allavginstanfreq));
allavgamp=allavgamp(~isnan(allavgamp));


end


function val=getdose(str)

if contains(str,'baseline','IgnoreCase',true)
    val=0;
elseif contains(str,'DAMGO','IgnoreCase',true)
    index1=strfind(str,'DAMGO');
    index2=strfind(str,'nM');
    str=char(str);
    val=str((index1+6):(index2-1));
    val=str2num(val);
else
    error("incorrect data labels")
end

end

function [avg_instan_freq,avg_freq,avg_instan_amp]=getrecdata(datastruct,dur)

%get burst times, amplitudes
bursttime=datastruct.rawBurstTimes;
burstamp=datastruct.rawBurstAmplitudes;
instanfreq=1./diff([0 bursttime]);
starttime=datastruct.filtTime(1); %in sec
endtime=datastruct.filtTime(end); %in sec

%get index for data only in the last 5 minutes
datalen=dur;
newstarttime=endtime-dur;
if newstarttime<starttime
    newstarttime=starttime;
    datalen=endtime-newstarttime;
end
dataindex=find(bursttime>newstarttime);

%avg data from last 5 minutes of recording
if isempty(dataindex)
    avg_instan_freq=0;
    avg_freq=0;
    avg_instan_amp=0;
else
    %get instan freq metrics
    avg_instan_freq=mean(instanfreq(dataindex));
    avg_freq=length(dataindex)/datalen;
    %get amplitude metrics
    avg_instan_amp=mean(burstamp(dataindex));
end

end

