function [wCoefs,peakidxs,peakInfo,removedPeaks]=MassSpecWavelet(data,time,Fs,Fsmax,SNR_Th,wavelet,analysiscond)
%this function is based on MassSpecWavelet code written by and based on
%"Improved peak detection in mass spectrum by incorporating continuous 
%wavelet transform-based pattern matching "
%(Du et al 2006);

%%ie numsamps=scale*16+1
maxscale=round(Fs/Fsmax);
scales=1:1:maxscale; %think about this...
wCoefs=cwt(data,scales,wavelet);

if strcmp(analysiscond,"manual")
figure;
g=gcf;
xlabel('scales');
title('Continuous Waveform Coefficients');

image(wCoefs,'CDataMapping','scaled');
colorbar

[X,Y]=ginput(1);

[~,~,th] = isoutlier(newdata,'mean','ThresholdFactor', 3); %defines the 
%threshold from the mean for outliers, ie peaks
[~,peakidxs]=findpeaks(newdata,'MinPeakProminence',th);

elseif strcmp(analysiscond,"automated")
    
    %get local maxima for coefficients at adjacent scales    
    localMax = getLocalMaximumCWT_2(wCoefs,scales);
    
    %identify ridge lines by connecting local maxima at adjacent scales
    [ridgeList,removedPeaks]= getRidge(localMax,scales);
    
    %extract peaks from ridgelines, and get peak attributes such as peakSNR
    [peakidxs,peakInfo,removedPeaks]=identifyMajorPeaks(data,time,ridgeList,removedPeaks,wCoefs,scales,SNR_Th);
        
end
    

end









