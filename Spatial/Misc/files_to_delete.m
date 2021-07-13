function files_to_delete(ei_C)
ind = 0;
for ii = 1:length(ei_C)
    ei = ei_C{ii};
    for pp = 1:length(ei.plane)
        thispFolder = ei.plane{pp}.folder;
        files = dir(thispFolder);
        for jj = 1:length(files)
            if ~isempty(strfind(files(jj).name,'0_motion'))
                ind = ind + 1;
                disp(ind);
                disp(files(jj).name);
                file_name = fullfile(thispFolder,files(jj).name);
                delete(file_name);
            end
        end
        n = 0;
    end
end
disp(ind)