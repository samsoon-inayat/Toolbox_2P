function mSigO = fixMSig(xs,mSig)
mSigP = [0 mSig 0];
for ii = 1:length(mSig)
    if isnan(mSig(ii))
        for jj = (ii+1):length(mSig)
            if ~isnan(mSig(jj))
                break;
            end
        end
        if ii == 1
            txs = [(ii-1) jj];
            txs1 = (ii-1):jj;
            temp = interp1(txs+1,mSigP(txs+1),txs1+1);
            mSig(txs1(2:end)) = temp(2:end);
            continue;
        end
        txs = [(ii-1) jj];
        txs1 = (ii-1):jj;
        temp = interp1(txs,mSig(txs),txs1);
        mSig(txs1) = temp;
    end
end
mSigO = mSig;