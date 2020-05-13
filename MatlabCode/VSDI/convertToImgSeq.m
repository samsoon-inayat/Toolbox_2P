function dfbf = convertToImgSeq(dfbf1,mask)

dfbf = zeros(size(mask,1),size(mask,2),size(dfbf1,2));
dfbf = reshape(dfbf,size(mask,1)*size(mask,2),size(dfbf1,2));
maskI = find(mask);
dfbf(maskI,:) = dfbf1;
dfbf = reshape(dfbf,size(mask,1),size(mask,2),size(dfbf1,2));
dfbf = applyMask(dfbf,mask);