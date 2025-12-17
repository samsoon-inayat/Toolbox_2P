function animal = get_exp_info(mD,animal_list,date_list)

rdata_dir = mD.rdata_dir;
pdata_dir = mD.pdata_dir;
adata_dir = mD.adata_dir;

face_cam = 'face';
pupil_cam = 'pupi';
paws_cam = 'video';
mat_file = 'recording';

for an = 1:length(animal_list)
    tanimal = animal_list{an};
    tdate = date_list{an};
    animal(an).ID = tanimal;
    animal(an).date = tdate;
    animal(an).rdir = fullfile(fullfile(rdata_dir,tanimal),date_list{an});
    animal(an).pdir = fullfile(fullfile(pdata_dir,tanimal),date_list{an});
    if ~exist(animal(an).pdir,'dir')
        mkdir(animal(an).pdir);
    end
    animal(an).adir = fullfile(fullfile(adata_dir,tanimal),date_list{an});
    if ~exist(animal(an).adir,'dir')
        mkdir(animal(an).adir);
    end
    files = dir(fullfile(animal(an).rdir, '*'));

    % Initialize fields
    animal(an).video.h264.face  = '';
    animal(an).video.h264.pupil = '';
    animal(an).video.h264.paws  = '';
    animal(an).mat              = '';
    
    for fn = 1:length(files)
        fname = files(fn).name;
        fullfname = fullfile(animal(an).rdir,fname);
        % Skip folders
        if files(fn).isdir
            continue;
        end
    
        % ---- H264 CAMERAS ----
        if endsWith(fname, '.h264', 'IgnoreCase', true)
            if startsWith(fname, 'face', 'IgnoreCase', true)
                animal(an).video.h264.face = fullfname;
            elseif startsWith(fname, 'pupi', 'IgnoreCase', true)
                animal(an).video.h264.pupil = fullfname;
            elseif startsWith(fname, 'video', 'IgnoreCase', true)
                animal(an).video.h264.paws = fullfname;
            end
        end
    
        % ---- MAT FILE ----
        if endsWith(fname, '.mat', 'IgnoreCase', true)
            animal(an).mat = fullfname;
        end
    end

end

