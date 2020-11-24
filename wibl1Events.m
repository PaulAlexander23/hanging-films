function [value,isTerminal,direction] = wibl1Events(t,y,timerID,timeout)

    [value1, isTerminal1, direction1] = eTimeout(timerID,timeout);
    [value2, isTerminal2, direction2] = eNan(y);
    [value3, isTerminal3, direction3] = eDewetWIBL1(y);

    value = [value1, value2, value3];
    isTerminal = [isTerminal1, isTerminal2, isTerminal3];
    direction = [direction1, direction2, direction3];
end

