function sig = normalizeSignal(sig,dim)

if isvector(sig)
    minSig = min(sig);
    sig = sig - minSig;
    maxSig = max(sig);
    sig = sig/maxSig;
    return;
end

if ~exist('dim','var')
    dim = 1;
end

if dim == 1
    minSig = min(sig);
    sig = sig-repmat(minSig,size(sig,1),1);
    maxSig = max(sig);
    sig = sig./repmat(maxSig,size(sig,1),1);
    return;
end

if dim == 2
    minSig = min(sig,[],2);
    sig = sig-repmat(minSig,1,size(sig,2));
    maxSig = max(sig,[],2);
    sig = sig./repmat(maxSig,1,size(sig,2));
    if sum(isnan(sig(:))) > 0
        sig(isnan(sig)) = 0;
    end
    return;
end