function out = lap(H,L)
    dxx = diff_ps(H, 2, L(1));
    dyy = diff_ps(H', 2, L(2))';
    out = dxx + dyy;
end