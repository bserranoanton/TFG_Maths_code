%This code is desgined to simulate system 5.1.
%By Belen Serrano Anton
%Created 03/03/2020
%Last Modified 31/03/2020

syms t_cell(t) p(t) 

%Constants
a = 3.5;
b = 0.29;

k = 0.3;
lambda = 0.9;

t0 = 0; 
tf = 20; 
dt_cell=diff(t_cell,t);

%Initial Conditions
c1 = 3;     %P(0)
c2 = 1;     %T(0)
c3 = 0;     %T'(0)
y0 = [c1 c2 c3];
eq1 = diff(t_cell,t,2) == -k*t_cell + lambda*p;
eq2 = diff(p,t) == a*p - b*t_cell*p;

vars = [t_cell(t); p(t)];
[V,S] = odeToVectorField([eq1,eq2]);

M = matlabFunction(V,'vars', {'t','Y'});
interval = [t0 tf];  %Time interval    
%Impose a nonnegativity constraint 
option1 = odeset('Events', @myEvent);
option2 = odeset('NonNegative',2); %T >= 0
%opt = [option1 option2];

%ySol = ode45(M,interval,y0, option1);
[t,ySol,te,ye,ie] = ode45(M,interval,y0, option1);
tValues = linspace(interval(1),interval(2),1000);

%--------------------PATOGENO-----------------------------------------------
%yValuesP = deval(ySol,tValues,1); %number 1 denotes first position: pathogen

yValuesP = ySol(:,1);

%Hacemos que una vez que el patógeno es 0 no vuelva a reproducirse
flag_defeated = yValuesP < 0.01;
yValuesP(flag_defeated) = 0;

min_p_cell = min(yValuesP);
index_time_min_p_cell = yValuesP == min_p_cell;
time_min_p_cell = min(t(index_time_min_p_cell));

%Hacemos que una vez que el patógeno es 0 no vuelva a reproducirse
% if(min_p_cell < 0.01)
%     yValuesP(min(find(index_time_min_p_cell==1)):end) = 0;
% end

max_p_cell = max(yValuesP);
index_time_max_p_cell = yValuesP == max_p_cell;
time_max_p_cell = t(index_time_max_p_cell);

%--------------------CELULAS T----------------------------------------------
%yValuesT = deval(ySol,tValues,2);
yValuesT = ySol(:,2);
min_t_cell = min(yValuesT);
index_time_min_t_cell = yValuesT == min_t_cell;
tValues(1) = 500; %Evita que el tiempo donde las células T son 0 sea 0
time_min_t_cell = min(t(index_time_min_t_cell));
tValues(1) = 0;

%yValuesT(find(index_time_min_t_cell==1):end) = 0;
 
max_t_cell = max(yValuesT);
index_time_max_t_cell = yValuesT == max_t_cell;
time_max_t_cell = t(index_time_max_t_cell);


%Plot results
figure 
hold on
[hA2] = plot(t,yValuesP,'r','LineWidth', 1);  %Pathogen
 
hold on
[hA1] = plot(t,yValuesT,'b','LineWidth', 1);  %T cells

hold on
%plot max
[hM1] = plot(time_max_t_cell,max_t_cell,'o','MarkerFaceColor','b', 'MarkerEdgeColor', 'b');
hold on
[hM2] = plot(time_max_p_cell,max_p_cell,'o','MarkerFaceColor','r', 'MarkerEdgeColor', 'r');
%plot min
hold on
[hM3] = plot(time_min_t_cell,min_t_cell,'s','MarkerFaceColor','b', 'MarkerEdgeColor', 'b');
hold on
[hM4] = plot(time_min_p_cell,min_p_cell,'s','MarkerFaceColor','r', 'MarkerEdgeColor', 'r');


disp(max_p_cell);
disp(time_max_p_cell);
disp(time_min_p_cell);

disp('----');

disp(max_t_cell);
disp(time_max_t_cell);
disp(time_min_t_cell);

xlabel('Tiempo');  ylabel('Número de células');
legend([hA2,hA1,hM1,hM2,hM3,hM4],'Patógeno','Células T', 'Max celulas T', ...
   'Max celulas patogeno', 'Min celulas T', 'Min celulas patogeno');

