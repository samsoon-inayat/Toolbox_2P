function ImgSeq = applyTemporalFilter(ImgSeq)

[width,height,nFrames] = size(ImgSeq);

Gauss1d = LPF_Gauss_25Hz; Lhalf = round((length(Gauss1d)-1)/2);
parfor r = 1:width
    for c = 1:height
        sig = ImgSeq(r,c,:); sig = sig(:);
        ImgSeq(r,c,:) = conv(sig,Gauss1d,'same');
    end
end

ImgSeq = ImgSeq(:,:,Lhalf+1:nFrames-Lhalf);