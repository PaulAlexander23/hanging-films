function [value, isTerminal, direction] = eTimeout(timerID, timeout)
    value = toc(timerID) - timeout;
    isTerminal = 1;
    direction = 0;
end
