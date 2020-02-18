function [localMax]=getLocalMaximumCWT(wCoefs,scales)
%this function is based on getlocalMaximumCWT code written by and based on
%"Improved peak detection in mass spectrum by incorporating continuous 
%wavelet transform-based pattern matching "
%(Du et al 2006);
%goal is to:
% Identify the local maximum of each column in 2-D CWT coefficients matrix by
% using a slide window. The size of slide window linearly changes from the 
% coarse scale (bigger window size) to detail scale. The scale of CWT 
% increases with the column index

sz=size(wCoefs); nrows=sz(1); ncols=sz(2);

localMax=zeros(nrows,ncols);
for i=1:ncols
    scale=scales(i);
    tempdata=wCoefs(:,i);
    [~,~,th] = isoutlier(tempdata,'mean','ThresholdFactor', 2); %defines the threshold from the mean for outliers, ie peaks
    [pks,locs]=findpeaks(tempdata,'MinPeakProminence',th);
    localMax(locs,i)=1;
end

end

