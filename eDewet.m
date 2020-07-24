function [value, isTerminal, direction] = eDewet(y)
    value = min(y,[],[1,2]);
    isTerminal = 1;
    direction = 0;
end
