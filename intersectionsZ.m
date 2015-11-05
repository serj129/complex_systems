function [x0,y0] = intersectionsZ(x1,y1,x2,y2)
x1 = x1(:);
y1 = y1(:);
x2 = x2(:);
y2 = y2(:);
x0 = [];
y0 = [];
for i = 1 : length(x1)-1
  a1 = [x1(i);y1(i)];
  a2 = [x1(i+1);y1(i+1)]; 
  for j = 1: length(x2)-1
    b1 = [x2(j);y2(j)];
    b2 = [x2(j+1);y2(j+1)]; 
    [xx,yy] = intersectLines(a1,a2,b1,b2);
    x0 = [x0,xx];
    y0 = [y0,yy];
  end
end
