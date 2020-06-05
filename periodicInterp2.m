function Vq = periodicInterp2(X, Y, V, Xq, Yq, method)
    %Assuming equal spacing between points x0, x1 and x2
    if nargin < 6, method = 'linear'; end

    if isvector(X)
        X = [X(2)-2*X(1), reshape(X, 1, [])];
        Y = [Y(2)-2*Y(1), reshape(Y, 1, [])];
    else
        X = [X(:,2)-2*X(:,1), X];
        X = [X(end,:); X];
        Y = [Y(2,:)-2*Y(1,:); Y];
        Y = [Y(:,end), Y];
    end
    V = [V(:,end), V];
    V = [V(end,:); V];

    Vq = interp2(X, Y, V, Xq, Yq, method);
end
