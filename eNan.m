function [value, isTerminal, direction] = eDewetWIBL1(y)
    value = double(~any(isnan(y),[1,2]));
    isTerminal = 1;
    direction = 0;
end

