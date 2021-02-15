%This code is desgined to simulate figure 5.2.
%By Belen Serrano Anton
%Created 03/03/2020
%Last Modified 31/03/2020

function res = macro_nond_toler_into(a_star, b_star)
syms t_cell(t) p(t) 


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

% Impose a nonnegativity constraint 
option2 = odeset('NonNegative',2); %T >= 0

ySol = ode45(M,interval,y0, option2);
tValues = linspace(interval(1),interval(2),1000);
yValuesP = deval(ySol,tValues,1); 

%If pathogen molecules are least than 0.01 we consider that pathogen has
%been totally defeated
if(min(yValuesP) <= 0.01)
    res = 1; %Intolerance
else
    res = 0; %Tolerance
end



