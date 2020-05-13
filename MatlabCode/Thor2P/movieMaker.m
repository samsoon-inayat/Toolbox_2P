function movieMaker(ei,varargin)
if nargin == 2
    frames = varargin{1};
end
if nargin == 3
    frames = varargin{1};
    type = varargin{2};
end

if ~exist('frames')
    frames = 1:100;
end
if ~exist('type')
    type = 'avi';
end

    
if strcmp(type,'avi')
rows = ei.pixelY;
cols = ei.pixelX;
hWaitBar = waitbar(0,sprintf('Processed Frames -'));
fid = fopen(ei.mcRawFile,'r');

V = VideoWriter(ei.movieFile,'Uncompressed AVI');
open(V);
for jj = 1:length(frames)
    waitbar(jj/length(frames),hWaitBar,sprintf('Processing Frame %d',jj));
    fseek(fid, (frames(jj)-1)*rows*cols*2, 'bof');
    I1 = fread(fid,rows*cols,'uint16=>double',0,'l');
    I1 = reshape(I1,rows,cols)';
    I1 = I1/max(I1(:));
    I1 = imadjust(I1);
    RGB = insertText(I1,[10 10],sprintf('%d',frames(jj)),'BoxColor','black','TextColor','white') ;
    writeVideo(V,RGB);
  end
fclose(fid);
close(V);
close(hWaitBar) ;
end

if strcmp(type,'tif')
    if exist(ei.movieFileTif)
        delete(ei.movieFileTif);
    end
    rows = ei.pixelY;
    cols = ei.pixelX;
    hWaitBar = waitbar(0,sprintf('Processed Frames -'));
    fid = fopen(ei.mcRawFile,'r');
    for jj = 1:length(frames)
        waitbar(jj/length(frames),hWaitBar,sprintf('Processing Frame %d',jj));
        fseek(fid, (frames(jj)-1)*rows*cols*2, 'bof');
        I1 = fread(fid,rows*cols,'uint16=>double',0,'l');
        I1 = reshape(I1,rows,cols)';
        I1 = I1/max(I1(:));
%         I1 = imadjust(I1);
%         RGB = insertText(I1,[10 10],sprintf('%d',frames(jj)),'BoxColor','black','TextColor','white') ;
        try
                if jj > 1
                    imwrite(I1,ei.movieFileTif,'WriteMode','append');
                else
                    imwrite(I1,ei.movieFileTif);
                end
        catch
            display('An error occured but I am trying again');
            if jj > 1
                    imwrite(I1,ei.movieFileTif,'WriteMode','append');
                else
                    imwrite(I1,ei.movieFileTif);
                end
        end
    
      end
    fclose(fid);
    close(hWaitBar) ;
end