function err = calculateError(x, t, data1, data2)
    domain = Domain(x); 

    M = domain.shape(1);
    N = domain.shape(2);

    M1 = data1.solution.domain.shape(1);
    M2 = data2.solution.domain.shape(1);
    nos = size(data1.solution.y,1)/M1; % Number of surfaces
    err = nan(M*nos, N, length(t));
    y1 = zeros(M*nos, N);
    y2 = zeros(M*nos, N);

    for n = 1:length(t)

        y1t = permute( ...
            interp1(data1.solution.t, permute(data1.solution.y, [3,1,2]), t(n)), ...
            [2,3,1]);
        
        for j = 1:nos
            y1(1+(j-1)*M:j*M,:) =  domain.interp(data1.solution.domain.x, y1t(1+(j-1)*M1:j*M1,:));
        end

        y2t = permute( ...
            interp1(data2.solution.t, permute(data2.solution.y, [3,1,2]), t(n)), ...
            [2,3,1]);

        for j = 1:nos
            y2(1+(j-1)*M:j*M,:) =  domain.interp(data2.solution.domain.x, y2t(1+(j-1)*M2:j*M2,:));
        end

        err(:,:,n) = y1 - y2;
    end
end
