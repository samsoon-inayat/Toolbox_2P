function temp_iqw
filename = 'question.xlsx';

[num,strings,raw] = xlsread(filename,1,'A2:B109');

n = 0;

[r,p] = corrcoef(num)

figure(100);
scatter(num(:,1),num(:,2));
xlabel('1st Column');ylabel('2nd Column');
format_axes_b(gca);
set(gca,'FontSize',12);
title(sprintf('Mutiple R = %.4f, p-value = %3.4e',r(1,2),p(1,2)));