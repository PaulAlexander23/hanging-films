function dy = compute_diff_ps(y, degree, xL, suppression)
    %DIFF_PSEUDO_SPECTRAL Uses the pseudo-spectral method to differentiate
    %   Detailed explanation goes here
    if nargin < 4
       suppression = 1e-13;
    end

    % Transform into fourier space
    yF = fft(y);

    N = size(yF,1)/2;

    % Determine k in matlab form
    %k = fftshift([0,-N+1:N-1])';
    k = [0:N-1, 0, 1-N:-1]' * 2*pi/xL;

    % Prior suppression
    yF(abs(yF)<suppression) = 0;

    % Apply pseudo-spectral differentiation
    dyF = (1i*k).^degree.*yF;

    % Posterior suppression
    % dyF(abs(dyF) < suppression*N*2) = 0 ;
    % dyF(abs(dyF) < suppression*max(abs(dyF))) = 0 ;

    % Transform back into real space
    dy = ifft(dyF);

end
