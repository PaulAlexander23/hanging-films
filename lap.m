function out = lap(H,L)
    dxx = compute_diff_ps(H, 2, L(1));
    dyy = compute_diff_ps(H', 2, L(2))';
    out = dxx + dyy;
end