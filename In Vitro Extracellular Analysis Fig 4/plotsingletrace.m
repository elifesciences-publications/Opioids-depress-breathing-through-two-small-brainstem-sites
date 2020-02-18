function plotsingletrace(datafile,ch,fs)
% the purpose of this code is to take in extracellular data, and plot
% integrated traces of voltage over time

%% initialize parameters for filter and plotting

% set filter for bandpass
lowerbound=100;
upperbound=4000;
nyquist=fs/2;
[B,A]=butter(2,[lowerbound/nyquist upperbound/nyquist],'bandpass'); %2nd order bandpass butterworth filter
initpath=pwd;

%% load in all data from file, and convert these data into time and voltage signal

[filepath,name,ext]=fileparts(datafile);

%get rec condition name, this changes based on your naming scheme
temp=strsplit(name,'_'); date=[temp{2} '/' temp{3} '/' temp{1}]; pclampid=temp{4};
recCond=strjoin(temp(5:end),'_');

switch ext
    case '.txt'
        
        fileID = fopen(datafile,'r');
        formatSpec = '%f';
        A=fscanf(fileID,formatSpec);
        
        DataVar.(recCond).raw=A;
        
        time=A(1:2:end)'; Vtrace=A(2:2:end)';
        Vtrace=Vtrace-mean(Vtrace); %normalize
        
        % remove first 10 seconds, bc pclamp generates artifacts at
        % beginning of recording
        time=time(10*fs:end);
        Vtrace=Vtrace(10*fs:end);
        
        %this may be case specific for our recordings, bc sometimes they
        %are recorded in either mV or V
        if max(Vtrace)<0.05
            Vtrace=Vtrace*1000;
        end
        
        rawmin=min(Vtrace)-0.1;
        rawmax=max(Vtrace)+0.1;
        
    case '.mat'
        
        DataVar.(recCond).raw=load(datafile);
        
        time=DataVar.(recCond).raw.c001_Time;
        Vtrace=DataVar.(recCond).raw.(ch);
        Vtrace=Vtrace-mean(Vtrace); %normalize
        
        % remove first 10 seconds, bc pclamp generates artifacts at
        % beginning of recording
        time=time(10*fs:end);
        Vtrace=Vtrace(10*fs:end);
        
        %this may be case specific for our recordings, bc sometimes they
        %are recorded in either mV or V
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

intmin=min(intTrace)-0.5;
intmax=max(intTrace)+0.5;

%% plot integrated,filtered data

figure;
g=gcf;
plot(newTime,intTrace,'color','k');
title(strcat(strrep(name,'_',' '),': Integrated Trace'));
xlabel('sec');ylabel('mV');
ylim([intmin intmax]);
hold on;

%% plot raw data

figure;
f=gcf;
plot(time,Vtrace,'color','k');
title(strcat(strrep(name,'_',' '),': Raw Trace'));
xlabel('sec');ylabel('mV');
ylim([rawmin rawmax]);
hold on;

%% save figures

rootfilepath=erase(filepath,'\Raw Data');
cd(rootfilepath)

if ~exist(strcat(rootfilepath,'\','Figures'),'dir')
    mkdir('Figures')
end
cd('Figures')

%go to channel folder
figpath=pwd;
if ~exist(strcat(figpath,'\',ch),'dir')
    mkdir(ch)
end
cd(ch)

savefig(g,strcat('integrated_',name));
saveas(g,strcat('integrated_',name,'.pdf'));
close(g)
savefig(f,strcat('raw_',name));
saveas(f,strcat('raw_',name,'.pdf'));
close(f)

cd(initpath) %go back to original path

end
