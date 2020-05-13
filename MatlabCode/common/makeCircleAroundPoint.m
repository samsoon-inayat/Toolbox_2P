function [xs,ys] = makeCircleAroundPoint (center,radius)

xo = center(1); yo = center(2);
sx = xo-radius+1; ex = xo+radius+1;
sy = yo-radius+1; ey = yo+radius+1;
[X,Y] = meshgrid(sx:ex,sy:ey);
ang = 0:pi/100:2*pi;
circum = [(radius*cos(ang') + xo) (radius*sin(ang') + yo)];
xs = circum(:,1);
ys = circum(:,2);

