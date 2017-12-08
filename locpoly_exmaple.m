clear;
x = -pi/2:0.01:pi/2;
y = sin(x);
[m0] = locpoly(y, x);
degree = 2; derivative = 2; width = 0.03; kernel ='box';
[m0, dm0, d2m0] = locpoly(y, x, [] ,degree, derivative, width, kernel);
plot(x,y,'*')
hold on 
plot(x, m0, 'm');
hold off
% estimate 1-order derivative, has good performance
plot(x, cos(x), '*')
hold on
plot(x, dm0, 'm') 
hold off
% estimate 2-order derivative, has bad effect
plot(x, -sin(x), '*')
hold on 
plot(x, d2m0, 'm')
hold off