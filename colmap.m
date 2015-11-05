clear all; clc;clf;
global tau1 = 1;
global tau2 = 1;
global nu_1 = 0.009;
global nu_2 = 0.01;
global beta0_1 = 0.042;
global beta0_2 = 0.042;
global q0 = 1;
global q1 = 1.7;
global C = 200;
global xi1 = 0.0006;
global xi2 = 0.0006;


%S-shape for first system N^(1)_2 from \beta^(0)_1
figure(1)
axes('linewidth',2,'fontsize',12, 'box', 'off','fontname', 'Arial');
xlabel('$\beta^{(1)}_0$');
ylabel('$\nu^{(1)}$');
hold on; 
beta0_1_min = 0;
beta0_1_max = 0.14;
beta0_1 = beta0_1_min:0.01:beta0_1_max;
nu_1_min = 0.0084;
nu_1_max = 0.012;
nu_1 = nu_1_min:0.0001:nu_1_max;
k= 0.000001:0.01:5;

for l = 1 : length(nu_1)
for j = 1 : length(beta0_1)
for i = 1 : length(k)

capitalB1(i) =(k(i)^3*tau2*q1/tau1/q0+k(i)^2*(1-nu_1(l)*C*tau2*q1)+k(i))/(tau1*q0+k(i)*tau1*q0+k(i)^2*tau2*q1);
N2overC1(i) =k(i)^2*tau2*q1/(tau1*q0+k(i)*tau1*q0+k(i)^2*tau2*q1);
N22(i)=(capitalB1(i)-beta0_1(j))/xi1;

capitalB2(i) =(k(i)^3*tau2*q1/tau1/q0+k(i)^2*(1-nu_2*C*tau2*q1)+k(i))/(tau1*q0+k(i)*tau1*q0+k(i)^2*tau2*q1);
N2overC2(i) =k(i)^2*tau2*q1/(tau1*q0+k(i)*tau1*q0+k(i)^2*tau2*q1);
N12(i)=(capitalB2(i)-beta0_2)/xi2;
end

[x0,y0,iout,jout]=intersections(N22,N2overC2*C,N2overC1*C,N12);
colmap(l,j)=length(y0);
%xx = ones(size(y0))*beta0_1(j);
%plot (xx,y0,'<');
end
end
imagesc([beta0_1_min beta0_1_max],[nu_1_min nu_1_max], colmap);
colorbar
print -depslatexstandalone '-S600,450' 'colormap.tex'
