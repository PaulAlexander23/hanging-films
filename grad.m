function out = grad(H, L)
    dx = compute_diff_ps(H, 1, L(1));
    dy = compute_diff_ps(H', 1, L(2))';
    out = zeros([size(H),2]);
    out(:,:,1) = dx;
    out(:,:,2) = dy;
end