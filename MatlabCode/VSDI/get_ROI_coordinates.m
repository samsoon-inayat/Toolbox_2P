function out = get_ROI_coordinates(ags,animal_group,animalNumber)

folder = ags(animal_group).animals{animalNumber}.root_folder_surjeet;

fileName = makeName('ROI_coordinates.mat',folder);

temp = load(fileName);
out = temp.ROI_coordinates;


