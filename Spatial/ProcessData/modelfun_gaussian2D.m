function z = modelfun_gaussian2D(b,xy)

dc = b(1);
A = b(2);
mux = b(3);
sigx = b(4);
muy = b(5);
sigy = b(6);


A1 = b(7);
mux1 = b(8);
sigx1 = b(9);
muy1 = b(10);
sigy1 = b(11);

A2 = b(12);
mux2 = b(13);
sigx2 = b(14);
muy2 = b(15);
sigy2 = b(16);


x = xy(:,1);
y = xy(:,2);

z = dc + A * exp(-(((x-mux)/sigx).^2 + ((y-muy)/sigy).^2)) + A1 * exp(-(((x-mux1)/sigx1).^2 + ((y-muy1)/sigy1).^2)) + A2 * exp(-(((x-mux2)/sigx2).^2 + ((y-muy2)/sigy2).^2));