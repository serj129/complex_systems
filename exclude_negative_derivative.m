function [x0,y0] = exclude_negative_derivative(x1,y1)

xd = diff(x1);
xx = x1(2:length(x1));
yy = y1(2:length(x1));
%x0 = xx(find(xd >= 0));
%y0 = yy(find(xd >= 0));
yy(find(xd < 0)) = NaN;
x0 = xx;
y0 = yy;
end