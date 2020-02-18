
function [inspPeak, expPeak, inspDur, expDur, inspVt, expVt, breathVt]=getbreathvals(breathmat)
%get critical features of individual breaths such as inspiratory peak
%(inspPeak), expiratory peak (expPeak), inspiratory duration (inspDur),
%expiratory duration (expDur), inspiratory tidal volume (inspVt),
%expiratory tidal volume (expVt), and total breath volume (breathVt)

%inputs: cell array (breathMat) of 1xn breaths, each containing array of
%voltage values of variable length for each breath

numbreaths=length(breathmat);

%initialize output arrays
inspVt=nan(numbreaths,1);
expVt=nan(numbreaths,1);
breathVt=nan(numbreaths,1);
inspDur=nan(numbreaths,1);
expDur=nan(numbreaths,1);
inspPeak=nan(numbreaths,1);
expPeak=nan(numbreaths,1);

for numbreath=1:numbreaths
    
    breath=breathmat{numbreath};
    breathtime=0.001:0.001:length(breath)/1000;
    zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0); % Returns Zero-Crossing Indices Of Argument Vector
    zx=zci(breath);
    if isempty(zx)
        continue
    end
    
    %get breath durations
    id=zx(1);%samples/ms
    inspDur(numbreath,1)=id;
    ed=length(breath)-zx(1);%samples/ms
    expDur(numbreath,1)=ed; 
    
    %get breath tidal volumes and total volume
    ivt= trapz(breathtime(1:zx(1)),breath(1:zx(1))); %in ml
    inspVt(numbreath,1)= ivt; 
    evt= trapz(breathtime(zx(1):end),breath(zx(1):end)); %in ml
    expVt(numbreath,1)= evt; 
    bvt=ivt+evt; %net ml
    breathVt(numbreath,1)=bvt;
    
    %get peak vals    
    ip=min(breath);
    inspPeak(numbreath,1)=ip; 
    ep=max(breath);
    expPeak(numbreath,1)=ep; 
    
end

end



