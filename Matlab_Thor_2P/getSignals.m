function sig = getSignals(ei,sigType)
if strcmp(lower(sigType),'caraw')
    for ii = 1:length(ei.areCells)
        sig(ii,:) = ei.signals(ii,:);
    end
end

if strcmp(lower(sigType),'ca')
    for ii = 1:length(ei.areCells)
        sig(ii,:) = ei.deconv.caSigAll{ii};
    end
end

if strcmp(lower(sigType),'sp')
    for ii = 1:length(ei.areCells)
        sig(ii,:) = ei.deconv.spSigAll{ii};
    end
end