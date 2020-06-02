function bidishift = get_bidi_shift(ei,planeNumber,nFrames)

if ei.zFastEnable
    ii = planeNumber;
    frameNumsA = ii:(ei.zSteps+1):(5000*(ei.zSteps+1));
else
    frameNumsA = 1:5000;
end

frameNumS = 1:nFrames:5000;
frameNumE = nFrames:nFrames:5000;
ind = 1;
frameNums = frameNumsA(frameNumS(ind):frameNumE(ind));
frames = getFramesFromRaw(ei,frameNums);

% aFrame = mean(frames{1},3);
aFrame = max(frames,[],3);
[rows,cols] = size(aFrame);
rows1 = floor(rows/3);
cols1 = floor(cols/3);
xlimsB = [1 cols];
ylimsB = [1 rows];
ylimsZ = [rows1 (2*rows1)];
xlimsZ = [cols1 (2*cols1)];
bidishift = 0;
xlimsC = xlimsB; ylimsC = ylimsB;
while 1
    disp('Press ''a'' for average image, ''m'' for max image, ''s'' for saving');
    disp('Press ''f'' for getting more frames, ''z'' for toggling zoom');
    disp('Press ''q'' for generating error and quitting');
    hf = figure(1000);clf;
    set(hf,'Position',get(0,'Screensize')); 
    subplot 121
    imagesc(aFrame);colormap gray;axis equal;
    xlim(xlimsC); ylim(ylimsC);
    title(sprintf('BiDiPhase = %d',bidishift));
    saFrame = ShiftBiDi(bidishift,aFrame,size(aFrame,1),size(aFrame,2));
    subplot 122;
    imagesc(saFrame);colormap gray;axis equal;
    xlim(xlimsC); ylim(ylimsC);
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
    if keyVal == 122 % z
        if isequal(xlimsC,xlimsB)
            xlimsC = xlimsZ; ylimsC = ylimsZ;
        else
            xlimsC = xlimsB; ylimsC = ylimsB;
        end
    end
    if keyVal == 113 % q
        error
        break;
    end
    if keyVal == 97 % a
        aFrame = mean(frames,3);
    end
    if keyVal == 109 % m
        aFrame = max(frames,[],3);
    end
    if keyVal == 102 % f
        ind = ind + 1;
        frameNums = frameNumsA(frameNumS(ind):frameNumE(ind));
        tframes = getFramesFromRaw(ei,frameNums);
        frames = cat(3,frames,tframes);
        aFrame = mean(frames,3);
    end
end
close(hf);


