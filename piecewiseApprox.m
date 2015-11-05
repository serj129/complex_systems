function [x0,y0] = piecewiseApprox(x,y)
%x = 0:.01:4;
%y=5*cos(10.*x)+x.^3-2*x.^2-6*x+10;
deriv_y = diff(y)./diff(x);
xd = x(2:length(x));
yd = y(2:length(y));
delay_mult = deriv_y(1:length(deriv_y)-1).*deriv_y(2:length(deriv_y));
critical = xd(find(delay_mult < 0));
x0=[x(1),critical,x(length(x))];
y0=[y(1),yd(find(delay_mult < 0)),y(length(y))];

%y0 = fitPiecewiseLinearFunction(x, y, x0);


%plot(x,y)
%hold
%plot(xd(find(delay_mult < 0)),y(find(delay_mult <0)),'s');
%plot(x0,p,'o');
