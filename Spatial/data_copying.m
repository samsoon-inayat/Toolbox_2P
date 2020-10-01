% ei15_C = loadContextsResponses(ei15_C,[1 1],[1 1 1]);
% ei10_C = loadContextsResponses(ei10_C,[1 1],[1 1 1]);


old_folder = 'Z:\homes\brendan.mcallister\2P\Processed_Data\OLD';

filesToCopy = {'contexts.mat','contexts_trials_markers.mat'};

TT = ET15_C;
for ii = 1:size(TT,1)
    animal_id = TT{ii,1};
    pyProcessedFolder = TT{ii,7};
    ind = findstr(pyProcessedFolder{1},num2str(animal_id));
    slash_ind = findstr(pyProcessedFolder{1}(1:ind),'\');
    animal_folder = pyProcessedFolder{1}((slash_ind(end)+1):end);
    for jj = 1:length(filesToCopy)
        old_file_name = fullfile(old_folder,animal_folder,'suite2p\plane0\post_suite2p_matlab',filesToCopy{jj});
        new_file_name = fullfile(pyProcessedFolder{1},'suite2p\plane0\post_suite2p_matlab',filesToCopy{jj});
        if exist(old_file_name,'file')
            disp(old_file_name)
            copyfile(old_file_name,new_file_name);
        else
            disp(sprintf('Non-existent - %s',old_file_name));
        end
        old_file_name = fullfile(old_folder,animal_folder,'suite2p\plane1\post_suite2p_matlab',filesToCopy{jj});
        new_file_name = fullfile(pyProcessedFolder{1},'suite2p\plane1\post_suite2p_matlab',filesToCopy{jj});
        if exist(old_file_name,'file')
            disp(old_file_name)
            copyfile(old_file_name,new_file_name);
        else
            disp(sprintf('Non-existent - %s',old_file_name));
        end
    end
end
