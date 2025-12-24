function animal = load_dlc_labeled_filenames(animal)
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

            % 2. Define the directory to search
            searchDir = animal(an).pdir;
            
            % 3. Construct the search pattern (BaseName + labeled)
            % Use * as a wildcard to allow for text before, between, or after the strings
            pattern = fullfile(searchDir, ['*' base '*labeled*']);
            
            % 4. Find the files
            matchingFiles = dir(pattern);

            if isempty(matchingFiles)
                fprintf('Not Found file: for %s\n', cam);
                continue;
            end

            if length(matchingFiles) > 1
                fprintf('Multiple files: for %s\n', cam);
                continue;
            end
            
            % 5. Access the names
            for i = 1:length(matchingFiles)
                fprintf('Found file: %s\n', matchingFiles(i).name);
            end


            out_file = fullfile(animal(an).pdir, matchingFiles(i).name);
            animal(an).video.mp4_labeled.(cam) = out_file;
        end
    end
end