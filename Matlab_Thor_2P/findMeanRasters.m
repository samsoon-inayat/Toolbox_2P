function meanR = findMeanRasters(R,trials)

if ~exist('trials','var')
    for cn = 1:size(R,3)
        meanR(cn,:) = applyGaussFilt(nanmean(R(:,:,cn)),5);
    end
else
    if length(trials) == 1
        for cn = 1:size(R,3)
            meanR(cn,:) = applyGaussFilt(R(trials,:,cn),5);
        end
    else
        for cn = 1:size(R,3)
            meanR(cn,:) = applyGaussFilt(nanmean(R(trials,:,cn)),5);
        end
    end
end