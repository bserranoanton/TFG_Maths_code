%When T = 0 we stop the simulation
function [position, isterminal, direction] = myEvent(t,ySol)
    position   = ySol(2);
    isterminal = 1;   % Stop the integration
    direction  = 0;
end