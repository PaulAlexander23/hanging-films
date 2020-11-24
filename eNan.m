function [value, isTerminal, direction] = eDewetWIBL1(y)
    value = double(~any(any(isnan(y))));
    isTerminal = 1;
    direction = 0;
end

