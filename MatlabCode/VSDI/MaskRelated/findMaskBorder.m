function [x,y] = findMaskBorder(Mask)
% Finding Mask border and order the point.
% The function is written based on the angles order.

border = edge(Mask,'sobel',0.1);
[y,x] = find(border == 1);

cx = mean(x);
cy = mean(y);

a = atan2(y-cy, x-cx);
[~, order] = sort(a);
x = x(order); x(end+1) = x(1);
y = y(order); y(end+1) = y(1);