function save_mean_img(aei)

for ii = 1:length(aei)
    ei = aei{ii};
    for pl = 1:length(ei.plane)
        [ii pl]
        [ops,~,~] = get_ops(ei,pl);
        mimg = ops.meanImg - min(ops.meanImg(:));
        mimg = mimg/max(mimg(:));
        folder = ei.plane{pl}.folder;
        fileName = fullfile(folder,'mean_image_cellpose.jpg');
        fileName = sprintf('mimg_cp_%d_%d.jpg',ii,pl);
        imwrite(mimg,fileName);
    end
end
return;
