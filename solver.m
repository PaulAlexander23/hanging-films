function [y, t] = solver(pdefun, t, x, y0, method, timestepper)
    %SOLVER Computes the numerical solution up to tFinal
    %   Detailed explanation goes here
    
    shape = size(y0);
    y0 = reshape(y0, [prod(shape),1]);
    
    function [F,J] = func(pdefun,t,x,y,method,shape)
        if nargout == 1
            F = reshape(...
                pdefun(t, x, reshape(y,[shape,numel(y)/prod(shape)]), method),...
                [prod(shape),numel(y)/prod(shape)]);
        elseif nargout == 2
            [f,j] = pdefun(t, x, reshape(y,[shape,numel(y)/prod(shape)]), method);
            F = reshape(f, [prod(shape),numel(y)/prod(shape)]);
            J = j;
        end
    end
    
    y = timestepper(@(t,y) func(pdefun,t,x,y,method,shape),t,y0);

    if isstruct(y)
        t = y.x;
        y = y.y;
    end
    y = squeeze(reshape(y,[shape,length(t)]));
end