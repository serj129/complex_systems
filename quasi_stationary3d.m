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

%n12 = [0.01,0.03,0.1,0.3,1,3,5,10,20,30,50,100];
Tmin = 0;
Tmax = 200;
dT = 0.01;
T = Tmin:dT:Tmax;
F_0 = [100;99;1;66;66;67];
%{      
figure(1)
hold on;
for beta = -0.05 : 0.003 : 0.2
        beta0_1=beta;
        F = lsode('model', F_0 , T);
        plot(beta0_1,F(length(F),3),'o');
        F_0 = F(length(F),:);
end

F_0 = [50;50;100;66;66;67];
for beta = 0.2 : -0.003 : -0.05
   % for i = 1 : length(n12)
        beta0_1=beta;
        %beta0_2=beta;
        %F_0 = [(C-n12(i))/2;(C-n12(i))/2;n12(i);180;18;2];
        %F_0 = [(C-n12(i)-50);50;n12(i);75;75;50];
        F = lsode('model', F_0 , T);
        plot(beta0_1,F(length(F),3),'x');
        F_0 = F(length(F),:);
    %end
end

figure(2)
hold on;
beta0_1 = 0.042;
F_0 = [50;50;100;66;66;67];
for nu = 0.03 : -0.0005 : 0.001
        nu_1 = nu;
        F = lsode('model', F_0 , T);
        plot(nu_1,F(length(F),3),'o');
        F_0 = F(length(F),:);
end

F_0 = [100;99;1;66;66;67];
for nu = 0.001 : 0.0005 : 0.03
        nu_1 = nu;
        F = lsode('model', F_0 , T);
        plot(nu_1,F(length(F),3),'x');
        F_0 = F(length(F),:);
end
%}
%lsode 3d diagram
%{

figure(3)
hold on;
xlabel('$\beta^{(1)}_0$');
ylabel('$\nu^{(1)}_0$');
zlabel('$N_2^{(1)}$');
F_0 = [100;99;1;66;66;67];

for nu = 0.001 : 0.0006 : 0.02
  nu_1 = nu;
  nu_2 = nu;
  for beta = -0.05 : 0.006 : 0.2       
    beta0_1 = beta; 
    beta0_2 = beta;
        F = lsode('model', F_0 , T);
        plot3(beta0_1,nu_1,F(length(F),3),'.','linewidth',5);
        F_0 = F(length(F),:);
  end
end


%}

%analytical diagram
%S-shape for first system N^(1)_2 from \beta^(0)_1

figure(5)
axes('linewidth',2,'fontsize',12, 'box', 'off','fontname', 'Arial');
xlabel('$\nu^{(1)}_0$');
zlabel('$\beta^{(1)}$');
ylabel('$N_2^{(1)}$');
hold on;

beta0_1 = -0.05 : 0.003 : 0.2;
nu = 0.001 : 0.0003 : 0.02;
k= 0.001:0.005:1.5;

for kk = 1 : length(nu)
for j = 1 : length(beta0_1)
for i = 1 : length(k)
capitalB1(i) =(k(i)^3*tau2*q1/tau1/q0+k(i)^2*(1-nu(kk)*C*tau2*q1)+k(i))/(tau1*q0+k(i)*tau1*q0+k(i)^2*tau2*q1);
N2overC1(i) =k(i)^2*tau2*q1/(tau1*q0+k(i)*tau1*q0+k(i)^2*tau2*q1);
N22(i)=(capitalB1(i)-beta0_1(j))/xi1;
capitalB2(i) =(k(i)^3*tau2*q1/tau1/q0+k(i)^2*(1-nu(kk)*C*tau2*q1)+k(i))/(tau1*q0+k(i)*tau1*q0+k(i)^2*tau2*q1);
N2overC2(i) =k(i)^2*tau2*q1/(tau1*q0+k(i)*tau1*q0+k(i)^2*tau2*q1);
N12(i)=(capitalB2(i)-beta0_2)/xi2;
end

[N22_,N2overC2C]=exclude_negative_derivative(N22,N2overC2*C);
[N12_,N2overC1C]=exclude_negative_derivative(N12,N2overC1*C);
%[xS0,yS0]=intersections(N22_,N2overC2C,N2overC1C,N12_);
[xS0,yS0]=intersections(N22,N2overC2*C,N2overC1*C,N12);
xSx = ones(size(yS0))*beta0_1(j);
nuSnu = ones(size(yS0))*nu(kk);
plot3(nuSnu,yS0,xSx,'.','markersize',5,'linewidth',5);
%plot3(nuSnu,xSx,yS0,'linewidth',2);
end
end

figure(6)
ty = linspace (-8, 8, 41)';
tx =  linspace (-8, 8, 41)';
[xx, yy] = meshgrid (tx, ty);
r = sqrt (xx .^ 2 + yy .^ 2) + eps;
tz = sin (r) ./ r;
mesh (tx, ty, tz);
