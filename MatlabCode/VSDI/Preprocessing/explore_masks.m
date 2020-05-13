function explore_masks(ags,pags)

for gg = 1:length(ags)
    thisGroup = ags(gg);
    pthisGroup = pags(gg);
    for aa = 1:length(thisGroup.animals)
        thisAnimal = thisGroup.animals{aa};
        pthisAnimal = pthisGroup.animals{aa};
        winopen(thisAnimal.root_folder);
        winopen(pthisAnimal.root_folder);
        n = 0;
    end
end
