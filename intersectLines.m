function [x0,y0] = intersectLines(a1,a2,b1,b2)
x0 =[];
y0 =[];
 
r = a2-a1;
s = b2-b1;

crossProduct = @(a,b) (a(1)*b(2)-a(2)*b(1));

uNumenator = feval(crossProduct,b1-a1,r);
denominator = feval(crossProduct,r,s);
if (uNumenator == 0 && denominator == 0)
  % there are cases with overlapping or touching
  return
end
if (denominator == 0)
  return
end

u = uNumenator/denominator;
t = feval(crossProduct,b1-a1,s)/denominator;

if (u>=0 && u<=1 && t>=0 && t<=1)
  point = a1+t*r;
  x0 = point(1);
  y0 = point(2);
  return
end

