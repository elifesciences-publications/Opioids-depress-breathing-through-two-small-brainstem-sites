function [subTime,sumTrace]=movtrapz(window,time,Trace)
%compute integral for each window (in samples)
index0=1;
sumTrace=[];
subTime=[];
while index0<=(length(time)-window)
    index1=index0+window;
    miniTrace=Trace(index0:index1);
    sumTrace=[sumTrace trapz(miniTrace)];
    subTime=[subTime time(index1)];
    index0=index1+1;
end
end

