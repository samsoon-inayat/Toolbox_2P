function colormap_ig

cm = colormap(gray);
cm = flipud(cm(1:size(cm,1),:));
colormap(cm);