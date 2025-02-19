function [mainMask,all_masks,mimgg] = make_resp_frame(ei,pl,allresp)
[ops,astat,arecells] = get_ops(ei,pl);
FS = 8;
micronsPerPixel = ei.thorExp.widthUM/ei.thorExp.pixelX;
mimg = ops.meanImgE;
mimgg = ops.meanImg - min(ops.meanImg(:));
mimgg = mimgg/max(mimgg(:));
maskZ = zeros(size(mimg));
mainMask = maskZ;
orcells = cell_list_op(allresp,[],'or',1);
for jj = 1:size(allresp,2)
  selCells = allresp{1,jj};
  allmask = maskZ;
  iii = 0;
  for ii = 1:length(astat)
%     [jj ii]
    if ~logical(arecells(ii,1)) % check if it is a cell
      continue;
    end
    iii = iii + 1;
    if ~selCells(iii) % check if cell is responsive
        
      continue;
    end
    stat = astat(ii);
    mask = maskZ;
    mask(stat{1}.ipix) = 1;
    mask = mask';
    allmask(mask==1) = 1;
  end
  all_masks(:,:,jj) = allmask;
  mainMask = mainMask + jj*allmask;
end

mainMask(mainMask>size(allresp,2)) = 0;