function T = reduce_table(T,animal_list,date_list)

inds = [];
for ii = 1:length(animal_list)
    for jj = 1:size(T,1)
        if strcmp(T{jj,1},animal_list{ii}) && strcmp(T{jj,2},date_list{ii})
            inds = [inds jj];
        end
    end
end
T = T(inds,:);