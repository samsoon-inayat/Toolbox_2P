function mainMask = com_resp_frame(ei,pl,allresp)
[ops,astat,arecells] = get_ops(ei,pl);
FS = 8;
micronsPerPixel = ei.thorExp.widthUM/ei.thorExp.pixelX;
mimg = ops.meanImgE;
maskZ = zeros(size(mimg));
mainMask = maskZ;
for jj = 1:size(allresp,2)
  for ii = 1:size(allresp,2)
    pop1 = allresp{1,jj};
    pop2 = allresp{1,ii};
    iii = 0;
    for ii = 1:length(astat)
        if ~logical(arecells(ii,1))
        continue;
      end
      iii = iii + 1;
      if ~selCells(iii)
        continue;
      end
      stat = astat(ii);
      stat{1}.x;
      mask = mask';
      allmask(mask==1) = 1;
    end
  end
  mainMask = mainMask + jj*allmask;
end