function animal = process_h264(animal,owr)
% PROCESS_H264  Convert .h264 camera videos to .mp4 and store in pdir.
%
%   animal = PROCESS_H264(animal)
%
%   For each animal(an), this function looks for:
%       animal(an).rdir                  % raw data dir
%       animal(an).pdir                  % processed data dir
%       animal(an).video.h264.face
%       animal(an).video.h264.pupil
%       animal(an).video.h264.paws
%
%   It converts each existing .h264 file to .mp4 (using ffmpeg),
%   stores the converted file in animal(an).pdir, and then records:
%       animal(an).video.mp4.face
%       animal(an).video.mp4.pupil
%       animal(an).video.mp4.paws
%
%   NOTE: Requires ffmpeg to be installed and on the system PATH.

if ~exist('animal','var')
    animal = evalin('base','animal')
end

if ~exist('owr','var')
    owr = 0;
end

    cams = {'face','pupil','paws'};

    for an = 1:numel(animal)

        % Skip if no h264 information
        if ~isfield(animal(an),'video') || ~isfield(animal(an).video,'h264')
            continue;
        end

        for c = 1:numel(cams)
            cam = cams{c};

            % If we don't have this camera, skip
            if ~isfield(animal(an).video.h264, cam)
                continue;
            end

            h264_name = animal(an).video.h264.(cam);

            if isempty(h264_name)
                continue;
            end

            % Full input path
            in_file = h264_name;
            if ~exist(in_file,'file')
                warning('process_h264:FileNotFound', ...
                    'Could not find %s for animal %d', in_file, an);
                continue;
            end

            % Output mp4 path (same base name)
            [~, base, ~] = fileparts(h264_name);
            mp4_name = [base '.mp4'];
            % mp4_name = [base '.avi'];
            out_file = fullfile(animal(an).pdir, mp4_name);
            % If mp4 already exists, skip conversion
            if exist(out_file, 'file') & owr == 0
                fprintf('Skipping (already exists): %s\n', mp4_name);
                animal(an).video.mp4.(cam) = out_file;
                continue
            end

            % Build ffmpeg command (copy stream, no re-encode)
            cmd = sprintf('ffmpeg -y -r 60 -i "%s" -c:v copy "%s"', in_file, out_file); %
            % cmd = sprintf(['ffmpeg -y -loglevel error -r 60 -i "%s" -vf "fps=60" -c:v libx264 -crf 10 -preset slow -pix_fmt yuv420p "%s"' ],in_file, out_file);
            
            % the above mentioned is for mp4

            % % the following is for avi because it turned out that mp4 files
            % % are slower with VideoReader
            % % cmd = sprintf('ffmpeg -y -r 60 -i "%s" -c:v mjpeg -q:v 3 copy "%s"', in_file, out_file);
            % cmd = sprintf(['ffmpeg -y -r 60 -i "%s" -vf "fps=60" -c:v mjpeg -q:v 1 -pix_fmt yuvj422p "%s"'], in_file, out_file);
            % % ffmpeg -i your_video.mp4 -c:v mjpeg -q:v 3 your_video_mjpeg.avi


            fprintf('Converting %s -> %s ... Please wait until I am done ... not showing progress of ffmpeg\r', in_file, out_file);
            status = system(cmd);

            if status ~= 0
                warning('process_h264:FFmpegError', ...
                    'ffmpeg failed for %s (animal %d, cam %s)', in_file, an, cam);
                continue;
            end

            % Store the mp4 file name (just the name, not the full path)
            animal(an).video.mp4.(cam) = mp4_name;
            % animal(an).video.avi.(cam) = mp4_name;
        end
    end
end