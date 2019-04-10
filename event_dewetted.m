function [value, isterminal, direction] = event_dewetted(~, y)
   
    value = 2 * any(y < 0,'all') - 1; % -1 if wetted, 1 if dewetted
    if any(isnan(y),'all'), value = 1; end % 1 if Nans in solution
    isterminal = 1; % Terminates
    direction = 0; % Allow any approach direction
    
end
