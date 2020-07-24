function [value, isTerminal, direction] = eDewet(y)
    value = min(min(y));
    isTerminal = 1;
    direction = 0;
end
