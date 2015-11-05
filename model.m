function dF = model(F,T)

%F1 = N0
%F2 = N1
%F3 = N2
%F4 = M0
%F5 = M1
%F6 = M2
global tau1;
global tau2;
global nu_1;
global nu_2;
global beta0_1;
global beta0_2;
global q0;
global q1;
global C;
global xi1;
global xi2;

dF = zeros(size(F));
dF(1) = -q0*(nu_1*F(3)+beta0_1+xi1*F(6))*F(1)+F(2)/tau1;
dF(2) =  q0*(nu_1*F(3)+beta0_1+xi1*F(6))*F(1)-F(2)/tau1-q1*(nu_1*F(3)+beta0_1+xi1*F(6))*F(2)+F(3)/tau2;
dF(3) =                                                 q1*(nu_1*F(3)+beta0_1+xi1*F(6))*F(2)-F(3)/tau2;
dF(4) = -q0*(nu_2*F(6)+beta0_2+xi2*F(3))*F(4)+F(5)/tau1;
dF(5) =  q0*(nu_2*F(6)+beta0_2+xi2*F(3))*F(4)-F(5)/tau1-q1*(nu_2*F(6)+beta0_2+xi2*F(3))*F(5)+F(6)/tau2;
dF(6) =                                               q1*(nu_2*F(6)+beta0_2+xi2*F(3))*F(5)-F(6)/tau2;


end
