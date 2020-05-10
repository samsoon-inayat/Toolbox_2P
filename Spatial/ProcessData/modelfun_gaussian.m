function y = modelfun_gaussian(b,x)

a1 = b(1);
b1 = b(2);
c1 = b(3);
% 
% a2 = b(4);
% b2 = b(5);
% c2 = b(6);
% 
% a3 = b(7);
% b3 = b(8);
% c3 = b(9);
% 
% 
% 
y =  a1*exp(-((x-b1)/c1).^2);% + a2*exp(-((x-b2)/c2).^2) + a3*exp(-((x-b3)/c3).^2);