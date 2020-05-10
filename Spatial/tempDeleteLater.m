function tempDeleteLater

ei = evalin('base','ei{1}');

cell_ids = ei.areCells';

spks1 = double(ei.tP.spks(cell_ids,:));

for ii = 1:length(cell_ids)
    spks2(ii,:) = ei.deconv.spSigAll{ii};
    
    figure(10);
    clf;
    plot(spks1(ii,:),'r');hold on;
    plot(spks2(ii,:),'b');
    pause
end

