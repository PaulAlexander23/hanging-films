function yhat = suppressHighWavenumbers(yhat, ratio)
    if nargin < 2, ratio = 0; end

    M = size(yhat, 1);
    N = size(yhat, 2);

    yhat(M/2 + (ceil(-M/2*ratio+1):floor(M/2*ratio)), :) = 0;
    yhat(:, N/2 + (ceil(-N/2*ratio+1):floor(N/2*ratio))) = 0;
end
