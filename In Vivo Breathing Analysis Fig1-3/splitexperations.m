function [ExpParams]=splitexperations(mouseid,AllBreaths,SaveCond)
%seperate expiratory period into 'expiration' and 'pause' based on a minimal flow threshold

%input: AllBreaths, struct containing fields of raw recording data in 
%different conditions, with AllBreaths.condition containing 1xn cells of breaths

fields=fieldnames(AllBreaths);

for numrec=1:length(fields)
    
    field=string(fields(numrec));
    breathmat=getfield(AllBreaths,field);
    
    %initialize arrays
    numbreaths=length(breathmat);
    Expdur=nan(numbreaths,1);
    Pause=nan(numbreaths,1);

    for numbreath=1:numbreaths
        
        breath=breathmat{numbreath};
        zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0); % Returns Zero-Crossing Indices Of Argument Vector
        zx=zci(breath);
        if isempty(zx)
            continue
        end
        
        %define pause as period after initiation of expiration where flow
        %drops below thresh
        thresh=0.5;
        exppeakind=find(breath==max(breath),1,'first');
        expend=find(breath(exppeakind:end)<=thresh,1,'first');
        if ~isempty(expend)
            expend=exppeakind+expend-1;
        else
            expend=length(breath);
        end
        
        Expdur(numbreath,1)=expend-zx(1);
        Pause(numbreath,1)=length(breath)-expend;

%%if you want to troubleshoot where it is defining pause, plot pauses 
%         subplot(1,2,2)
%         plot(breath,'k');
%         hold on; plot([expend expend],[-5 5],'--r')
%         
%         close(h)
        
    end
    
    ExpParams.(field).Expdur=Expdur;
    ExpParams.(field).Pause=Pause;

end

if SaveCond
    savestr=string(strcat(mouseid,'_ExpParams.mat'));
    save(savestr,'ExpParams');
end

end



