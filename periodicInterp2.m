function Vq = periodicInterp2(X, Y, V, Xq, Yq, method)
    %Assuming equal spacing between points x0, x1 and x2
    if nargin < 6, method = 'linear'; end

    V = [V(:,end), V];
    V = [V(end,:); V];

    if isvector(X)
        X = [X(2)-2*X(1), reshape(X, 1, [])];
        Y = [Y(2)-2*Y(1), reshape(Y, 1, [])];

        Vq = interp2(Y, X, V, Yq, Xq, method);
    else
        X = [X(:,2)-2*X(:,1), X];
        X = [X(end,:); X];
        Y = [Y(2,:)-2*Y(1,:); Y];
        Y = [Y(:,end), Y];

        Vq = interp2(X, Y, V, Xq, Yq, method);
    end
end
