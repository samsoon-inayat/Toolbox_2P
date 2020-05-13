function dfbyfo = PreProcessEV_HiLo_OnlyDfbyfo_N(hi,lo,no,Fs)
[width,height,~] = size(hi);
% Load no trials and use it to calc F0
F0 = no;
for r = 1:width
    parfor c = 1:height
        sig = F0(r,c,:);sig =sig(:);
        yhat = sig - locdetrend(sig,Fs,[0.3,0.1]);
        F0(r,c,:) = yhat;
    end
end
dfbyfo_hi = 100*(hi-F0)./F0; dfbyfo_hi(isnan(dfbyfo_hi))=0; dfbyfo_hi(isinf(dfbyfo_hi))=0;

dfbyfo_lo = 100*(lo-F0)./F0; dfbyfo_lo(isnan(dfbyfo_lo))=0; dfbyfo_lo(isinf(dfbyfo_lo))=0;

dfbyfo = (dfbyfo_hi + dfbyfo_lo)/2;

