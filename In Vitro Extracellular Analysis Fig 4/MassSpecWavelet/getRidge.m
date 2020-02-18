function [ridgeList,orphans]= getRidge(localMax,scales)
%this function is based on getlocalMaximumCWT code written by and based on
%"Improved peak detection in mass spectrum by incorporating continuous
%wavelet transform-based pattern matching "
%(Du et al 2006);
%goal is to:
% Return a list of ridge. As some ridges may end at the scale larger than 1,
% in order to keep the uniqueness of the ridge names, we combined the
% smallest scale of the ridge and m/z index of the peak at that scale
% together to name the ridges. For example the ridge name "1\_653" means
% the peak ridge ends at the CWT scale 1 with m/z index 653 at scale 1.

iInit=size(localMax);iInit=iInit(2);
stepsize=-1;
iFinal=1;
minWinSize=5;
gapTh=3; %how many times a peak can be 'not found' in a convolution before it is thrown out

maxInd_curr=find(localMax(:,iInit)>0); %get indices for largest scale peaks
nMz=size(localMax); nMz=nMz(1); %get num rows
nCol=size(localMax);nCol=nCol(2);

% identify all the peak pathes from the coarse level to the detailed level
%(higher scale values to lower scale values)
% only consider the shortest path (so its not winding around)
if nCol>1
    colInd=(iInit+stepsize):stepsize:iFinal; %everything but current index, iInit
else
    colInd=1;
end

%initialize struct
fieldnames={'peakNum' 'peakIdx' 'peakStatus' 'ridgeIdx'};
ridgeList=initializeStruct(fieldnames,maxInd_curr,true);

%orphan ridge is supposed to keep the ridges disconnected at certain scales
orphanRidgeList=[];
orphanRidgeName=[];
nlevel=length(colInd);

for j=1:nlevel %for every scaled convolution smaller than Iind (Iind is std of comparison)
    
    col=colInd(j);
    scale=scales(col);
    
    if length(maxInd_curr)==0 %bug?
        maxIndcurr=find(localMax(:,col)>0);
        error('bug in code')
    end
    
    winSize=scale*2+1;
    if winSize<minWinSize
        winSize=minWinSize;
    end
    
    %initialize new vals and remove vals before k iter
    selPeak=[];
    remove=[];
    
    %% find indices for where peaks line up in localmax(:,col), ie which
    %peaks are duplicated across scales
    
    for k=1:length(maxInd_curr) %for every 'peak' in Iind, determine if it is upheld by this scale
        ind=maxInd_curr(k); %index of peak in the Iind cond (ie largest scale)
        if ind-winSize<1 startval=1; else startval=ind-winSize; end
        if ind+winSize>nMz endval=nMz; else endval=ind+winSize; end
        ind_curr=find(localMax(startval:endval,col)>0) + startval - 1;
        %index of peak in window of next scale val - startval (so relative to peak in largest scale now)
        if length(ind_curr)==0
            ridgenum=find(cell2mat(cellfun(@(x) ismember(x,ind),{ridgeList.ridgeIdx},'UniformOutput',false))); %get index where peakIdx=idx
            status=ridgeList(ridgenum).peakStatus; %check current status
            if isempty(status) %no previous peaks found to match Iind peak here
                ridgeList(ridgenum).peakStatus=gapTh+1; % mark as a confirmed peak
                status=ridgeList(ridgenum).peakStatus; %check current status
            end
            if status>gapTh && scale>=2 % there is no peak at this val
                %mark this ridge for removal
                ridgenum=find(cell2mat(cellfun(@(x) ismember(x,ind),{ridgeList.ridgeIdx},'UniformOutput',false))); %get index where peakIdx=idx
                orphanRidgeList=[orphanRidgeList ridgeList(ridgenum).ridgeIdx]; %what the orig value was in the largest scale
                orphanRidgeName=[orphanRidgeName sprintf('%d_%d',col,ind)]; %col is which filter removed it, dupPeak is the duplicate value found
                remove=[remove ridgenum];
                continue
            elseif scale<2 %scale<2 ie scale=1 OR peak has been prev validated but not found now
                ind_curr=ind;
                ridgenum=find(cell2mat(cellfun(@(x) ismember(x,ind),{ridgeList.ridgeIdx},'UniformOutput',false))); %get index where peakIdx=idx
                ridgeList(ridgenum).peakStatus=status+1; % dock a point against it
            else %same but if status is not > gapTh, n
                ind_curr=ind;
                ridgenum=find(cell2mat(cellfun(@(x) ismember(x,ind),{ridgeList.ridgeIdx},'UniformOutput',false))); %get index where peakIdx=idx
                ridgeList(ridgenum).peakStatus=status+1; % dock a point against it
            end
        else %there is a value in this scale within the window size
            ridgenum=find(cell2mat(cellfun(@(x) ismember(x,ind),{ridgeList.ridgeIdx},'UniformOutput',false))); %get index where peakIdx=idx
            status=ridgeList(ridgenum).peakStatus; %check current status
            if isempty(status)
                ridgeList(ridgenum).peakStatus=0; % mark as a confirmed peak
            end
            if length(ind_curr)>=2 %choose which to use
                tempind_curr=abs(ind_curr-ind);
                ind_curr=ind_curr(find(tempind_curr==min(tempind_curr))); %take the value with the smallest absolute distance from ind
                ind_curr=ind_curr(1); %just in case they are the same dist away
            end
            %add this index value to the peakIdx for comparison
            ridgeList(ridgenum).peakIdx=[ridgeList(ridgenum).peakIdx ind_curr];
            %             ridgeList(ind).ridgeIdx=[ridgeList(ind).ridgeIdx ind_curr];
            if ~(length(ind_curr)==length(ridgenum))
                error('bug in code')
            end
            selPeak=[selPeak ind_curr]; %the ones we're adding
            
        end
    end
    
    %% update ridgeList
    
    % remove the disconnected lines from the current list
    if length(remove)>0
        ridgeList(remove)=[];
    end
    ridgeList=removeemptyindices(ridgeList);
    
    %resort struct by new ridgeidx
    ridgeList=resortbyridgeIdx(ridgeList);
    
    %renumber the peaks
    ridgeList=resetpeakNum(ridgeList);
    
    %update the ridgelist ids, based on selPeakIdx, which contains peaknum
    %val
    ridgeList=updateridgeIdx(ridgeList);
    
    %% add any additional peaks found in this new scale
    
    newPeak=[];
    currentPeaks=getridgeIdx(ridgeList);
    if j<nlevel/2 %halway through
        maxInd_next=find(localMax(:,col)>0); %get all peaks from this col
        %see which peaks are not already in selPeak ie recognized
        unSelPeak=setdiff(maxInd_next,selPeak)';
        %create new peak struct with these values
        if ~isempty(unSelPeak)
            newPeak=initializeStruct(fieldnames,unSelPeak,false);
            %update rigelist
            ridgeList=vertcat(ridgeList,newPeak);
            maxInd_curr=[currentPeaks unSelPeak];
        else
            maxInd_curr=[currentPeaks];
        end
    else
        maxInd_curr=[currentPeaks];
    end
    
    %% update ridgeList
    
    %resort struct by new ridgeidx
    ridgeList=resortbyridgeIdx(ridgeList);
    
    %renumber the peaks
    ridgeList=resetpeakNum(ridgeList);
    
    %% check for duplicated selected peaks and only keep the ones with the longest path
    
    currentPeaks=getridgeIdx(ridgeList);
    dupindices=findduplicates(currentPeaks);
    dupPeaks=unique(currentPeaks(dupindices)); %what is the value of the duplicates
    if length(dupPeaks)>0
        removeInd=[]; %there are duplicates
        for dupidx=1:length(dupPeaks) %in case there is multiple duplicates
            %find indices for each peak in selPeak
            dupPeak=dupPeaks(dupidx);
            %find where these values are in ridgeList
            ridgenum=find(cell2mat(cellfun(@(x) ismember(x,dupPeak),{ridgeList.ridgeIdx},'UniformOutput',false))); %get index where peakIdx=dupPeak
            
            %pick ridgenum with longest path
            [~,ridgenum]=getlongestPath(ridgeList,ridgenum);
            % peakidxnum=find(ridgeList(peaknum).peakIdx==dupPeak));
            
            %set these peaks for removal
            removeInd=[removeInd ridgenum];
            
            %add them to orphanRidgeList
            orphanRidgeList=[orphanRidgeList ridgeList(ridgenum).ridgeIdx]; %what the orig value was in the largest scale
            orphanRidgeName=[orphanRidgeName  {sprintf('%d_%d',col,dupPeak)}]; %col is which filter removed it, dupPeak is the duplicate value found
            
        end
        selPeak=unique(selPeak);
        ridgeList(removeInd)=[];
        ridgeList=removeemptyindices(ridgeList);
        
    end
    
    %% update ridgeList
    
    %resort struct by new ridgeidx
    ridgeList=resortbyridgeIdx(ridgeList);
    
    %renumber the peaks
    ridgeList=resetpeakNum(ridgeList);
    
    %% update orphaned peaks ie. invalidated peaks
    
    orphans.orphanRidgeList=orphanRidgeList;
    
end

end

function [longest,shorter]=getlongestPath(myStruct,indices)
temp=struct2cell(myStruct)';
temparr=temp(indices,2)';
temparr=cellfun(@(x) length(x),temparr);
longest=indices(find(temparr==max(temparr),1,'first'));
shorter=indices(indices~=longest);
end

function arr=getridgeIdx(myStruct)
temp=struct2cell(myStruct)';
arr=cell2mat(temp(:,4))';
end

function newStruct=initializeStruct(fieldnames,peakindices,indexcase)
numfields=length(fieldnames);
numindices=length(peakindices);

if numindices==0
    newStruct=struct(fieldnames{:});
    return
end

newStruct=cell(numindices,numfields);
if indexcase
    newStruct(:,1)=num2cell(1:numindices,numindices);
end
newStruct(:,2)=num2cell(peakindices);
newStruct(:,4)=num2cell(peakindices);
newStruct=cell2struct(newStruct,fieldnames,2);
end

function newStruct=updateridgeIdx(myStruct)
fields=fieldnames(myStruct);
temp=struct2cell(myStruct)';
subTemp=temp(:,2);
subTemp=num2cell(cellfun(@(x) x(end),subTemp));
temp(:,4)=subTemp;
newStruct=cell2struct(temp,fields,2);
end

function newStruct=removeemptyindices(myStruct)
fields=fieldnames(myStruct);
temp=struct2cell(myStruct)';
emptyidx=find(sum((cellfun('isempty',temp)),2)>0);
if ~isempty(emptyidx)
    temp(emptyidx,:)=[];
end
newStruct=cell2struct(temp,fields,2);
end

function newStruct=resetpeakNum(myStruct)
sz=length(myStruct);
fields=fieldnames(myStruct);
temp=struct2cell(myStruct)';
temp(:,1)=num2cell(1:1:sz);
newStruct=cell2struct(temp,fields,2);
end

function newStruct=resortbyridgeIdx(myStruct)
fields=fieldnames(myStruct);
temp=struct2cell(myStruct)';
oldOrder=cell2mat(temp(:,4)');
[newOrder,newIndices]=sort(oldOrder,'ascend');
if ~(isequal(oldOrder,newOrder))
    temp(:,1:4)=temp(newIndices,1:4);
end
newStruct=cell2struct(temp,fields,2);
end

function [indx]=findduplicates(x)
%find index of duplicate values in an array
[~,i,~] = unique(x,'first');
indx=find(not(ismember(1:numel(x),i)));
end
