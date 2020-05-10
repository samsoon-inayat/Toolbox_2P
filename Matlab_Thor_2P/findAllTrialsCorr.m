function [corrVal mCorrVal] = findAllTrialsCorr(signal)
corrVal = zeros(1,size(signal,1));
for ii = 1:size(signal,1)
    rr = 0;
    for jj = 1:size(signal,1)
%         if ii == jj
%             continue;
%         end
        sig1 = signal(ii,:);
        sig2 = signal(jj,:);
        corrVal(ii,jj) = rr + xcorr(sig1,sig2,0,'coeff');
    end
%     corrVal(ii) = rr/(size(signal,1)-1);
end

for ii = 1:size(signal,1)
    temp = corrVal(ii,:);
    temp(ii) = [];
    mCorrVal(ii) = mean(temp);
end