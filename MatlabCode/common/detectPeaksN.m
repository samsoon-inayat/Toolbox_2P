%Function to detect peaks in a signal
%Usage
%detectedPeakPositions = detectPeaks(signal);
%Input is a vector signal
%Ouput is nx2 matrix with the first column having indices in the vector
%where peaks exist and the second column values of those indices.

function detectedPeakPositions = detectPeaksN (signal)

signalInv = 1.01*max(signal) - signal;
detectedPeakPositions = detectPeaks(signalInv);
if ~isempty(detectedPeakPositions)
    detectedPeakPositions(:,2) = signal(detectedPeakPositions(:,1));
end

