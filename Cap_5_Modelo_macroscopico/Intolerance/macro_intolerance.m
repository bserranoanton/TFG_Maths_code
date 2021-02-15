%This code is desgined to simulate system 5.1.
%By Belen Serrano Anton
%Created 03/03/2020
%Last Modified 31/03/2020

syms t_cell(t) p(t) 

%Constants
a = 1.5;
b = 0.1;

k = 0.4;
lambda = 0.5;

t0 = 0; 
tf = 10; 
dt_cell=diff(t_cell,t);

%Initial Conditions
c1 = 3;     %P(0)
c2 = 0;     %T(0)
c3 = 0;     %T'(0)
y0 = [c1 c2 c3];
eq1 = diff(t_cell,t,2) == -k*t_cell + lambda*p;
eq2 = diff(p,t) == a*p - b*t_cell*p;

vars = [t_cell(t); p(t)];
[V,S] = odeToVectorField([eq1,eq2]);

M = matlabFunction(V,'vars', {'t','Y'});
interval = [t0 tf];  %Time interval    
%Impose a nonnegativity constraint 
option2 = odeset('NonNegative',2); %T >= 0



ySol = ode45(M,interval,y0, option2);
tValues = linspace(interval(1),interval(2),1000);
yValues = deval(ySol,tValues,1); 

%Plot results 
figure 
[hA2] = plot(tValues,yValues/max(yValues),'r','LineWidth', 1);  %Pathogen

hold on
yValues = deval(ySol,tValues,2); 
[hA1] = plot(tValues,yValues/max(yValues),'b','LineWidth', 1);  %T cells 
ylim([0,1]);

set(gca,'YTickLabel',[]); 
set(gca,'XTickLabel',[]);
xlabel('Tiempo');  ylabel('Numero de celulas');
legend([hA2,hA1],'Patogeno','Celulas T');

