%This code is desgined to simulate system 5.2.
%By Bel�n Serrano Ant�n
%Created 03/03/2020
%Last Modified 31/03/2020

syms t_cell(t) p(t) 

a_star = 1.1;
b_star = 0.01;

t0 = 0; 
tf = 9.5; 
dt_cell=diff(t_cell,t);

%Initial Conditions
c1 = 1;     %P(0)
c2 = 0;     %T(0)
c3 = 0;     %T'(0)
y0 = [c1 c2 c3];
eq1 = diff(t_cell,t,2) == -t_cell + p;
eq2 = diff(p,t) == a_star*p - b_star*t_cell*p;

vars = [t_cell(t); p(t)];
[V,S] = odeToVectorField([eq1,eq2]);


M = matlabFunction(V,'vars', {'t','Y'});
interval = [t0 tf];  %Time interval    

%Impose a nonnegativity constraint 
option2 = odeset('NonNegative',2); %T >= 0

ySol = ode45(M,interval,y0, option2);
tValues = linspace(interval(1),interval(2),1000);
yValues = deval(ySol,tValues,1); %number 1 denotes first position: pathogen

%Plot results
figure 
[hA2] = plot(tValues,yValues,'r','LineWidth', 1);  %Pathogen
 
hold on
yValues = deval(ySol,tValues,2); 
[hA1] = plot(tValues,yValues,'b','LineWidth', 1);  %T cells 

set(gca,'YTickLabel',[]); 
set(gca,'XTickLabel',[]);
xlabel('Tiempo');  ylabel('N�mero de c�lulas');
legend([hA2,hA1],'Pat�geno','C�lulas T');

