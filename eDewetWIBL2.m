function [value, isTerminal, direction] = eDewetWIBL2(y)
    value = min(min(y(1:end/3,:)));
    isTerminal = 1;
    direction = 0;
end

