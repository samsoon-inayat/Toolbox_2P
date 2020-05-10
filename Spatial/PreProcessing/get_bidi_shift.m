function bidishift = get_bidi_shift(ei,planeNumber,nFrames)

if ei.zFastEnable
    ii = planeNumber;
    getFrames = ii:(ei.zSteps+1):(nFrames*(ei.zSteps+1));
    frames{1} = getFramesFromRaw(ei,getFrames);
else
    frames{1} = getFramesFromRaw(ei,1:nFrames);
end

aFrame = mean(frames{1},3);
bidishift = 0;
lims = [100 500];
while 1
    hf = figure(1000);clf;
    set(hf,'Position',get(0,'Screensize')); 
    subplot 121
    imagesc(aFrame);colormap gray;axis equal;
    xlim(lims); ylim(lims);
    title(sprintf('BiDiPhase = %d',bidishift));
    saFrame = ShiftBiDi(bidishift,aFrame,size(aFrame,1),size(aFrame,2));
    subplot 122;
    imagesc(saFrame);colormap gray;axis equal;
    xlim(lims); ylim(lims);
    keyVal = getkey;%pause;
    if keyVal == 115 % s to return
%         if ei.zFastEnable
%             allbds(ii) = bidishift;
%         end
        break;
    end
    if keyVal == 29 % Up arrow
        bidishift = bidishift + 1;
    end
    if keyVal == 28 % down arrow
        bidishift = bidishift - 1;
    end
    if keyVal == 27 % down arrow
        pos = get(gcf,'Position');
        pos(3) = pos(3) + 10;
        pos(4) = pos(4) + 10;
        set(gcf,'Position',pos);
    end
    if keyVal == 108 % l
        lims = input('Enter limits ');
        xlim(lims); ylim(lims);
    end
    if keyVal == 113
        error
        break;
    end
end
close(hf);


