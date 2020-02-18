function [peakidxs,peakInfo,removedPeaks]=identifyMajorPeaks(data,time,ridgeList,removedPeaks,wCoefs,scales,SNRthresh,varargin)
%this function is based on MassSpecWavelet code (function identifyMajorPeaks)
% written by and based on
%"Improved peak detection in mass spectrum by incorporating continuous 
%wavelet transform-based pattern matching "
%(Du et al 2006);

nearbyWinSize=20; %1 sec
numPeaks=length(ridgeList);
removedPeaks.orphanMajorPeaks=[];

%attributes of all ridges
ridgeIdxArrays=getRidgeArrays(ridgeList); %list of all idx values for peak
ridgeLen=getRidgeLen(ridgeList);
peakStatus=getPeakStatus(ridgeList);
normPeakStatus=ridgeLen-peakStatus; %num confirmations- num dings
ridgeIdx=getRidgeIdx(ridgeList);
peakNum=1:numPeaks;
%remove all peaks that have low prob of being real (ie num confirmations <
%num dings)
removeIdx=find(peakStatus>3);
removeIdx=[removeIdx find(normPeakStatus<0)]; 
removeIdx=unique(removeIdx);
if ~isempty(removeIdx)
    removedPeaks.orphanMajorPeaks=[removedPeaks.orphanMajorPeaks ridgeIdx(removeIdx)];
    ridgeIdxArrays(removeIdx)=[];
    ridgeLen(removeIdx)=[];
    peakStatus(removeIdx)=[];
    normPeakStatus(removeIdx)=[];
    ridgeIdx(removeIdx)=[];
    numPeaks=length(ridgeLen);
    peakNum=1:numPeaks;
end


peakidxs=ridgeIdx;
peakvals=data(ridgeIdx);

%if 2 peaks are within nearbyWinSize of each other
peakdiff=diff(peakidxs);
tempremoveIdx=find(peakdiff<nearbyWinSize);
removeIdx=[];
if ~isempty(tempremoveIdx)
    for i=1:length(tempremoveIdx)
        indices=tempremoveIdx(i):(tempremoveIdx(i)+1);
        status=normPeakStatus(indices(1):indices(2)); % the two similar values
        statusvalidx=find(status==max(status)); %find which has smaller norm peak (more validated)
        statusvalidx=statusvalidx(1);
        index=indices(statusvalidx);
        removeIdx=[removeIdx index];
    end
end

%remove the problematic peaks
if ~isempty(removeIdx);
    removedPeaks.orphanMajorPeaks=[removedPeaks.orphanMajorPeaks ridgeIdx(removeIdx)];
    ridgeIdxArrays(removeIdx)=[];
    ridgeLen(removeIdx)=[];
    peakStatus(removeIdx)=[];
    normPeakStatus(removeIdx)=[];
    ridgeIdx(removeIdx)=[];
    peakidxs(removeIdx)=[];
    peakvals(removeIdx)=[];
    numPeaks=length(ridgeLen);
    peakNum=1:numPeaks;
end

%z-score wCoefs
newWCoefs=normalize_ib(wCoefs,2,'zscore');

%z-score data
newData=normalize_ib(data,1,'zscore');

%find peaks with low signal to remove
cfSNRval=[];
SNRval=[];
cfSNRval=mean(newWCoefs(ridgeIdx,:),2);
SNRval=newData(ridgeIdx);

removeIdx=[];
if ~isnan(SNRthresh)
removeIdx=find(SNRval<SNRthresh); %based on zscore from filtered data
end

%remove the problematic peaks
if ~isempty(removeIdx);
    removedPeaks.orphanMajorPeaks=[removedPeaks.orphanMajorPeaks ridgeIdx(removeIdx)];
    ridgeIdxArrays(removeIdx)=[];
    ridgeLen(removeIdx)=[];
    peakStatus(removeIdx)=[];
    normPeakStatus(removeIdx)=[];
    ridgeIdx(removeIdx)=[];
    cfSNRval(removeIdx)=[];
    SNRval(removeIdx)=[];
    peakidxs(removeIdx)=[];
    peakvals(removeIdx)=[];
    numPeaks=length(ridgeLen);
    peakNum=1:numPeaks;
end

peakInfo.peakNum=peakNum;
peakInfo.ridgeIdx=ridgeIdx;
peakInfo.peakStatus=peakStatus;
peakInfo.normPeakStatus=normPeakStatus;
peakInfo.ridgeLen=ridgeLen;
peakInfo.coeffSNR=cfSNRval;
peakInfo.SNR=SNRval;
peakInfo.peakIdx=peakidxs;
peakInfo.peakVal=peakvals;


end

function arr=getRidgeLen(myStruct)
temp=struct2cell(myStruct)';
arr=(temp(:,2))';
arr=cellfun(@(x) length(x),arr);
end

function arr=getPeakStatus(myStruct)
temp=struct2cell(myStruct)';
arr=cell2mat(temp(:,3))';
end


function arr=getRidgeIdx(myStruct)
temp=struct2cell(myStruct)';
arr=cell2mat(temp(:,4))';
end

function arr=getRidgeArrays(myStruct)
temp=struct2cell(myStruct)';
arr=(temp(:,2))';
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


