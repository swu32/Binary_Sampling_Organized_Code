% Raster Plot

function [] = plotRaster(spikeMat, tVec,color)
hold all;
spikeMat = spikeMat>0.5;% convert from array to logical index

for trialCount = 1:size(spikeMat,1)
%   if(~all(spikeMat(trialCount, :)))
    spikePos = tVec(spikeMat(trialCount, :));
    for spikeCount = 1:length(spikePos)
        plot([spikePos(spikeCount) spikePos(spikeCount)], ...
            [trialCount-0.4 trialCount+0.4], color);
    end
    
end
ylim([0 size(spikeMat, 1)+1]);