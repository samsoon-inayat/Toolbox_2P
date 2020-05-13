% x1 = -100:1:100; x2 = -100:1:100;
% [X1,X2] = meshgrid(x1, x2);
% sigma1 = 3;
% sigma2 = 3;
% scale1 = 3;
% scale2 = 6;
% sigma1 = scale1*sigma1;
% sigma2 = scale2*sigma2;
% Theta = 0;
% 
% a = ((cosd(Theta)^2) / (2*sigma1^2)) + ((sind(Theta)^2) / (2*sigma2^2));
% b = -((sind(2*Theta)) / (4*sigma1^2)) + ((sind(2*Theta)) / (4*sigma2^2));
% c = ((sind(Theta)^2) / (2*sigma1^2)) + ((cosd(Theta)^2) / (2*sigma2^2));
% 
% mu = [0 0];
% A = 1;
% F = A*exp(-(a*(X1 - mu(1)).^2 + 2*b*(X1 - mu(1)).*(X2 - mu(2)) + c*(X2 - mu(2)).^2));
% figure(1);clf;
% imagesc(F);


A = 1;
x0 = 0; y0 = 0;

sigma_x = 1;
sigma_y = 2;

[X, Y] = meshgrid(-5:.1:5, -5:.1:5);

for theta = 0:pi/100:pi
    a = cos(theta)^2/(2*sigma_x^2) + sin(theta)^2/(2*sigma_y^2);
    b = -sin(2*theta)/(4*sigma_x^2) + sin(2*theta)/(4*sigma_y^2);
    c = sin(theta)^2/(2*sigma_x^2) + cos(theta)^2/(2*sigma_y^2);

    Z = A*exp( - (a*(X-x0).^2 - 2*b*(X-x0).*(Y-y0) + c*(Y-y0).^2)) ;
    figure(1);clf;
    imagesc(Z)
end