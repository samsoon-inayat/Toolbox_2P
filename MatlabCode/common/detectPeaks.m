%Function to detect peaks in a signal
%Usage
%detectedPeakPositions = detectPeaks(signal);
%Input is a vector signal
%Ouput is nx2 matrix with the first column having indices in the vector
%where peaks exist and the second column values of those indices.

function detectedPeakPositions = detectPeaks (signal)

sizeSignal = max(size(signal));

peakPositionsIndex = 1;
if signal(2) < signal(1)
    detectedPeakPositions(peakPositionsIndex,1) = 1;
    detectedPeakPositions(peakPositionsIndex,2) = signal(1);
    peakPositionsIndex = peakPositionsIndex + 1;
end

movingIndex = 2;
while movingIndex < sizeSignal
    if signal(movingIndex-1) < signal(movingIndex) && signal(movingIndex + 1) < signal(movingIndex)
        detectedPeakPositions(peakPositionsIndex,1) = movingIndex;
        detectedPeakPositions(peakPositionsIndex,2) = signal(movingIndex);
        peakPositionsIndex = peakPositionsIndex + 1;
    end
    if signal(movingIndex-1) == signal(movingIndex) && signal(movingIndex) > 0
        plateauStart = movingIndex - 1;
        while movingIndex < sizeSignal
            if signal(movingIndex) ~= signal(movingIndex + 1) 
                break;
            end
            movingIndex = movingIndex + 1;
        end
        plateauEnd = movingIndex;
        if plateauStart > 1 && plateauEnd < sizeSignal
            if signal(plateauStart) > signal(plateauStart - 1) && signal(plateauEnd) > signal(plateauEnd + 1)
                totalPlateauValues = plateauEnd - plateauStart + 1;
                if rem(totalPlateauValues,2) == 0
                    plateauPeak = plateauStart + (totalPlateauValues/2);
                    detectedPeakPositions(peakPositionsIndex,1) = plateauPeak;
                    detectedPeakPositions(peakPositionsIndex,2) = signal(plateauPeak);
                    peakPositionsIndex = peakPositionsIndex + 1;
                    detectedPeakPositions(peakPositionsIndex,1) = plateauPeak-1;
                    detectedPeakPositions(peakPositionsIndex,2) = signal(plateauPeak-1);
                    peakPositionsIndex = peakPositionsIndex + 1;
                else
                    plateauPeak = plateauStart + floor(totalPlateauValues/2);
                    detectedPeakPositions(peakPositionsIndex,1) = plateauPeak;
                    detectedPeakPositions(peakPositionsIndex,2) = signal(plateauPeak);
                    peakPositionsIndex = peakPositionsIndex + 1;
                end
            end
        end
        if plateauStart == 1
            if signal(plateauEnd) > signal(plateauEnd + 1)
                totalPlateauValues = plateauEnd - plateauStart + 1;
                if rem(totalPlateauValues,2) == 0
                    plateauPeak = plateauStart + (totalPlateauValues/2);
                    detectedPeakPositions(peakPositionsIndex,1) = plateauPeak;
                    detectedPeakPositions(peakPositionsIndex,2) = signal(plateauPeak);
                    peakPositionsIndex = peakPositionsIndex + 1;
                    detectedPeakPositions(peakPositionsIndex,1) = plateauPeak-1;
                    detectedPeakPositions(peakPositionsIndex,2) = signal(plateauPeak-1);
                    peakPositionsIndex = peakPositionsIndex + 1;
                else
                    plateauPeak = plateauStart + floor(totalPlateauValues/2);
                    detectedPeakPositions(peakPositionsIndex,1) = plateauPeak;
                    detectedPeakPositions(peakPositionsIndex,2) = signal(plateauPeak);
                    peakPositionsIndex = peakPositionsIndex + 1;
                end
            end
        end
        if plateauEnd == sizeSignal
            if signal(plateauStart) > signal(plateauStart - 1)
                totalPlateauValues = plateauEnd - plateauStart + 1;
                if rem(totalPlateauValues,2) == 0
                    plateauPeak = plateauStart + (totalPlateauValues/2);
                    detectedPeakPositions(peakPositionsIndex,1) = plateauPeak;
                    detectedPeakPositions(peakPositionsIndex,2) = signal(plateauPeak);
                    peakPositionsIndex = peakPositionsIndex + 1;
                    detectedPeakPositions(peakPositionsIndex,1) = plateauPeak-1;
                    detectedPeakPositions(peakPositionsIndex,2) = signal(plateauPeak-1);
                    peakPositionsIndex = peakPositionsIndex + 1;
                else
                    plateauPeak = plateauStart + floor(totalPlateauValues/2);
                    detectedPeakPositions(peakPositionsIndex,1) = plateauPeak;
                    detectedPeakPositions(peakPositionsIndex,2) = signal(plateauPeak);
                    peakPositionsIndex = peakPositionsIndex + 1;
                end
            end
        end
        movingIndex = plateauEnd + 1;
    else
        movingIndex = movingIndex + 1;
    end
end

if signal(sizeSignal-1) < signal(sizeSignal)
    detectedPeakPositions(peakPositionsIndex,1) = sizeSignal;
    detectedPeakPositions(peakPositionsIndex,2) = signal(sizeSignal);
    peakPositionsIndex = peakPositionsIndex + 1;
end

if ~exist('detectedPeakPositions')
    detectedPeakPositions = [];
end