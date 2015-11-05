clear all; clc;clf;
global tau1 tau2 nu_1 nu_2 beta0_1 beta0_2 q0 q1 C xi1 xi2
tau1 = 1;
tau2 = 1;
nu_1 = 0.01;
nu_2 = 0.01;
beta0_1 = 0.042;
beta0_2 = 0.042;
q0 = 1;
q1 = 1.7;
C = 200;
xi1 = 0.001;
xi2 = 0.001;


k= 0.000001:0.01:1.2;

for i = 1 : length(k)
capitalB1(i) =(k(i)^3*tau2*q1/tau1/q0+k(i)^2*(1-nu_1*C*tau2*q1)+k(i))/(tau1*q0+k(i)*tau1*q0+k(i)^2*tau2*q1);
N2overC1(i) =k(i)^2*tau2*q1/(tau1*q0+k(i)*tau1*q0+k(i)^2*tau2*q1);
N22(i)=(capitalB1(i)-beta0_1)/xi1;
capitalB2(i) =(k(i)^3*tau2*q1/tau1/q0+k(i)^2*(1-nu_2*C*tau2*q1)+k(i))/(tau1*q0+k(i)*tau1*q0+k(i)^2*tau2*q1);
N2overC2(i) =k(i)^2*tau2*q1/(tau1*q0+k(i)*tau1*q0+k(i)^2*tau2*q1);
N12(i)=(capitalB2(i)-beta0_2)/xi2;
end

%{
[pwN2,pwBeta]=piecewiseApprox(N2overC2*C,N22);
[pwN21,pwBeta1]=piecewiseApprox(N2overC1*C,N12);

[x0,y0]=intersectionsZ(pwN2,pwBeta,pwN21,pwBeta1);
figure(5)
plot (x0,y0,'o','markersize',5,'linewidth',5);
figure(1)
axes('linewidth',2,'fontsize',12, 'box', 'off','fontname', 'Arial');
%axis([0,C,0,C],'equal');
xlabel('$N_2^{(2)}$');
ylabel('$N_2^{(1)}$');
hold on;
plot(N22,N2overC2*C,'linewidth',1);
plot(pwBeta,pwN2,'k','linewidth',2);
%plot(N22,N2overC2*C,'linewidth',3);
plot(N2overC1*C,N12,'linewidth',1,'color','k', 'linestyle','-.');
plot(pwN21,pwBeta1,'r','linewidth',2);
%}
figure(2)

axes('linewidth',2,'fontsize',12, 'box', 'off','fontname', 'Arial');
axis([0,150,0,150],'equal');
hold on;
ylabel('$x_1$');
xlabel('$y_1$');
%plot( capitalB1,N2overC1*C,'linewidth',5);
[xx,yy]=exclude_negative_derivative(N12,N2overC1*C);
%plot( N12,N2overC1*C,'linewidth',5);
plot( xx,yy,'linewidth',2);
plot( yy,xx,'linewidth',2);
nan = isnan(yy);
i1 = min(find(nan==1))-1;
i2 = max(find(nan==1))+1;
p1 = polyfit(yy(1:i1),xx(1:i1),2)
plot(polyval(p1,yy(1:i1)),yy(1:i1),':','linewidth',2)
p2 = polyfit(yy(i2:length(yy)),xx(i2:length(yy)),2)
plot(polyval(p2,yy(i2:length(yy))),yy(i2:length(yy)),':','linewidth',2)
print -depslatexstandalone '-S600,400' 'xy.tex'


%{
%S-shape for first system N^(1)_2 from \beta^(0)_1
figure(5)
axes('linewidth',2,'fontsize',12, 'box', 'off','fontname', 'Arial');
xlabel('$\beta^{(1)}_0$');
ylabel('$N_2^{(1)}$');
hold on;
beta0_1 = 0:0.001:0.1;
k= 0.000001:0.01:5;
for j = 1 : length(beta0_1)
for i = 1 : length(k)
capitalB1(i) =(k(i)^3*tau2*q1/tau1/q0+k(i)^2*(1-nu_1*C*tau2*q1)+k(i))/(tau1*q0+k(i)*tau1*q0+k(i)^2*tau2*q1);
N2overC1(i) =k(i)^2*tau2*q1/(tau1*q0+k(i)*tau1*q0+k(i)^2*tau2*q1);
N22(i)=(capitalB1(i)-beta0_1(j))/xi1;
capitalB2(i) =(k(i)^3*tau2*q1/tau1/q0+k(i)^2*(1-nu_2*C*tau2*q1)+k(i))/(tau1*q0+k(i)*tau1*q0+k(i)^2*tau2*q1);
N2overC2(i) =k(i)^2*tau2*q1/(tau1*q0+k(i)*tau1*q0+k(i)^2*tau2*q1);
N12(i)=(capitalB2(i)-beta0_2)/xi2;
end
[x0,y0]=intersections(N22,N2overC2*C,N2overC1*C,N12);
xx = ones(size(y0))*beta0_1(j);
plot (xx,y0,'.','markersize',5,'linewidth',5);

end

print -depslatexstandalone '-S600,450' 'n2beta0.tex'
%}
%{
figure(1)
axes('linewidth',2,'fontsize',12, 'box', 'off','fontname', 'Arial');
axis([0,0.25,0,1]);
xlabel('$B^{(1)}$');
ylabel('$N_2^{(1)}$');
hold on
plot(capitalB1,N2overC1);
%}
F_0 = [75;100;25;75;100;25];
Tmin = 0;
Tmax = 50;
dT = 0.01;
T = Tmin:dT:Tmax;
%{
F = lsode('model', F_0 , T);
figure(2);
hold on;
plot(T',F(:,1:3),'linewidth',3,'color','black');
figure(3);
hold on;
plot(T',F(:,4:6),'linewidth',3,'color','black');
%}
%{ phase portrait
f_0 = [
-10,10;
 10,-10;
 10,10;
 50,5;
 50,0;
  5,50;
  0,50;
 80,10;
 10,80;
 80,20;
 20,80;
 50,40;
 40,50;
 30,30;
 70,70;
 110,110
];
%{
figure(4)
axis([0,C,0,C],'equal');
hold on;
for i=1:length(f_0)
F_0 = [(C-f_0(i,1))/2,(C-f_0(i,1))/2,f_0(i,1),(C-f_0(i,2))/2,(C-f_0(i,2))/2,f_0(i,2)];
%F = lsode('model',F_0, T);
plot(F(:,6),F(:,3),'linewidth',2,'color','black');
endfor

plot(N22,N2overC2*C,'linewidth',1,'color','r', 'linestyle','-.');
plot(N2overC1*C,N12,'linewidth',1,'color','b', 'linestyle','-.');
%}
%}
%print -depslatexstandalone '-S600,450' 'n2overC.tex'
%{
figure(2)
axes('linewidth',2,'fontsize',12, 'box', 'off','fontname', 'Arial');
axis([0,5,0.001,0.7]);
xlabel('$\alpha_{20}C$');
ylabel('$N_2/C$');
hold on;
plot(alpha20C,n2oC,'linewidth',5,'color','k', 'linestyle','-');
plot(alpha20C1,n2oC,'linewidth',5,'color','k', 'linestyle','-');
plot(alpha20C2,n2oC,'linewidth',5,'color','k', 'linestyle','-');

text(1.6,0.2, sprintf('\\colorbox{white}{\\textbf{1}}'));
text(2.2,0.2, sprintf('\\colorbox{white}{\\textbf{2}}'));
text(2.5,0.2, sprintf('\\colorbox{white}{\\textbf{3}}'));
print -depslatexstandalone '-S600,450' 'n2fromalphaC.tex'

figure(3)
axes('linewidth',2,'fontsize',12, 'box', 'off','fontname', 'Arial');
axis([0,200,0,200]);%,'equal');
xlabel('$N^{(2)}_2$');
ylabel('$N^{(1)}_2$');
hold on;
plot((BETTA0.-0.001)./0.001,n2oC*300,'linewidth',5,'color','k', 'linestyle','-');
plot(n2oC*300,(BETTA0.-0.001)./0.001,'linewidth',5,'color','b', 'linestyle','--');

text(-1.2,0.65, sprintf('\\colorbox{white}{\\textbf{3}}'));
text(0.1,0.65, sprintf('\\colorbox{white}{\\textbf{1}}'));
text(-0.5,0.65, sprintf('\\colorbox{white}{\\textbf{2}}'));
print -depslatexstandalone '-S600,550' 'n2frombetta.tex'
%}
%{
figure('visible','off');
p=0.06:0.001:0.5;
k=7;
w=[0,0.5,1,2];
axes('linewidth',2,'fontsize',12, 'box', 'off','fontname', 'Arial');
axis([0,0.5,0,1]);
%plot(x,delta(x),'linewidth',2);
hold on;
for i=1:length(w)
plot (p,1.-1./(2*k.*p.-w(i)/2),'linewidth',3)
endfor

%plot([rho0], [0],'o','linewidth',10);
text(0.1, 0.15,'\colorbox{white}{$\rho_0$}', 'horizontalalignment','center');
%text(rho0, -0.2, \sprintf('\\colorbox{white}{$%.2f$}r', rho0), 'horizontalalignment','center');
title ("");  legend ("off");  grid();
xlabel('$p/r$');
ylabel('$i');
print -depslatex '-S450,300' 'ifromp.tex'
%}
