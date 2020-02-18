function [breathStartInd] = ek_breathSegmentationFx(breathArray, durThresh, inspAmpThresh, expAmpThresh)

% EK 06.25.2018
% segment breaths function
% input: breathArray (continuous breathing voltage trace), durThresh (breath duration threshold),
% inspAmpThresh (inspiration amplitude threshold), expAmpThresh (expiration amplitude threshold).

crossZero =[];
threshCnt = 0; % count up time from now until thresh - account for noise/minibreaths
isBreath = 0;
breathStartInd = 1;%[]; % idx of start (col1) and end (col2) of each breath (row)
crossZero = 1;
exp = 0;
for i = 2 : length(breathArray) % just find crossings first
    prev = sign(breathArray(i - 1));
    curr = sign(breathArray(i));
    if prev == -1 && curr == +1 % if trace crosses upwards past 0
        exp = 1;
    end
    if prev == +1 && curr == -1 % if trace crosses downward past 0 
            if exp % checks to make sure that an expiration happened between inspirations
                breathStartInd = [breathStartInd i]; 
                breathEndInd = breathStartInd(end) - 1; % also use start for setting breath end index
                exp = 0;
                if breathEndInd - breathStartInd(end - 1) > durThresh
                    if max(breathArray(breathStartInd(end - 1) : breathEndInd)) > expAmpThresh
                        if min(breathArray(breathStartInd(end - 1) : breathEndInd)) > inspAmpThresh %smaller than expected amplitude
                            breathStartInd = [breathStartInd(1 : end - 2) i];
                            exp = 1;
                        end
                    else
                        breathStartInd = [breathStartInd(1 : end - 2) i];
                        exp = 1;
                    end
                else
                    breathStartInd = [breathStartInd(1 : end - 2) i];
                    exp = 1;
                end
            end
    end
end

x = 1;




