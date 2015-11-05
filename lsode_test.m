clear all; clc;clf;
global tau1 tau2 nu_1 nu_2 beta0_1 beta0_2 q0 q1 C xi1 xi2;
tau1 = 1;
tau2 = 1;
nu_1 = 0.01;
nu_2 = 0.01;
beta0_1 = 0.042;
beta0_2 = 0.042;
q0 = 1;
q1 = 1.7;
C = 200;
xi1 = 0.0006;
xi2 = 0.0006;

n12 = [0.01,0.03,0.1,0.3,1,3,5,10,20,30,50,100];
Tmin = 0;
Tmax = 200;
dT = 0.01;
T = Tmin:dT:Tmax;
figure(1)
hold on;
for beta = -0.03 : 0.005 : 0.1
    for i = 1 : length(n12)
        beta0_1=beta;
        %F_0 = [(C-n12(i))/2;(C-n12(i))/2;n12(i);180;18;2];
        F_0 = [(C-n12(i)-50);50;n12(i);75;75;50];
        F = lsode('model', F_0 , T);
        plot(beta0_1,F(length(F),3),'o');
    end
end
%{
beta0_1 = -0.05:0.001:0.65;
k= 0.000001:0.01:1.5;
for j = 1 : length(beta0_1)
for i = 1 : length(k)
capitalB1(i) =(k(i)^3*tau2*q1/tau1/q0+k(i)^2*(1-nu_1*C*tau2*q1)+k(i))/(tau1*q0+k(i)*tau1*q0+k(i)^2*tau2*q1);
N2overC1(i) =k(i)^2*tau2*q1/(tau1*q0+k(i)*tau1*q0+k(i)^2*tau2*q1);
N22(i)=(capitalB1(i)-beta0_1(j))/xi1;
capitalB2(i) =(k(i)^3*tau2*q1/tau1/q0+k(i)^2*(1-nu_2*C*tau2*q1)+k(i))/(tau1*q0+k(i)*tau1*q0+k(i)^2*tau2*q1);
N2overC2(i) =k(i)^2*tau2*q1/(tau1*q0+k(i)*tau1*q0+k(i)^2*tau2*q1);
N12(i)=(capitalB2(i)-beta0_2)/xi2;
end


[xS0,yS0]=intersections(N22,N2overC2*C,N2overC1*C,N12);
xSx = ones(size(yS0))*beta0_1(j);
plot (xSx,yS0,'.','markersize',5,'linewidth',5);
end
%}