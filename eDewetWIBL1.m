function [value, isTerminal, direction] = eDewetWIBL1(y)
    value = min(y(1:end/2,:),[],[1,2]);
    isTerminal = 1;
    direction = 0;
end

