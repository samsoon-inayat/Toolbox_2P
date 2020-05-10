function colors = getColors(n,bk)

tempC = distinguishable_colors(n,bk);

for ii = 1:size(tempC,1)
    colors{ii} = tempC(ii,:);
end


