function yhat = suppressSmallModes(yhat, suppression)
    if nargin < 2, suppression = 1e-13; end

    yhat(abs(yhat)<suppression) = 0; % Absolute suppression
    %yhat(abs(yhat)<suppression * max(max(abs(yhat)))) = 0; Relative suppression
end
