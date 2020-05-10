add_to_path

%%
clear all
clc
[f,cName,allD] = getFolders;



%%
runSuite2P(f,'183628','2019-07-01','',1,1,NaN);

%%
animal_ind = find(strcmp(allD.animalList,'183628'));
winopen(allD.animal(animal_ind(1)).folder)

%% animal list

AD_Thy1_animals = [183224;183227;183228;183329];
Thy1_animals = [183633;183761;183745;183628;183762];

animalList = AD_Thy1_animals;
% animalList = Thy1_animals;
clc
for ii = 1:length(animalList)
    animal_ind = find(strcmp(allD.animalList,num2str(animalList(ii))));
    [allD.animalList(animal_ind)   allD.animal(animal_ind).dateList]
end