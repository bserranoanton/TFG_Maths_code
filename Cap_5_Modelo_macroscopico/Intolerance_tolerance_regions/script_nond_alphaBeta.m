%This code is desgined to simulate figure 5.2.
%By Belen Serrano Anton
%Created 03/03/2020
%Last Modified 31/03/2020

a = 0.1;
b = 0.1;

f1 = figure

xlabel('\beta^{*}');  ylabel('\alpha^{*}');
ylim([0,2.5]);
xlim([0,2.5]);

while (a <= 2.5)
    b = 0.1;
    while(b <= 2.5)
        %Result of system 5.2
        res =  macro_nond_toler_into(a, b);
        figure(f1)
        hold on
        if(res == 1) %Intolerance
            plot(b,a,'d','MarkerFaceColor','green', 'MarkerEdgeColor', 'green');
        else %Tolerance
            plot(b,a,'d','MarkerFaceColor','red', 'MarkerEdgeColor', 'red');
        end
        hold on
        b = b + 0.1;
    end
    a = a + 0.1;
end

