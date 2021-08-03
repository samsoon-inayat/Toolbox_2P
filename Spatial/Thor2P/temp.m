figure(1);
clf;
subplot 131
showCells(gca,ei15_A{3},1,[])
subplot 132
showCells(gca,ei15_AA{1},1,[])
subplot 133
showCells(gca,ei15_AA{2},1,[])

figure(1);clf;plot(0,0)
for ii = 1:length(ei15_C)
    showCells(gca,ei15_C{ii},1,[])
    pause;
end