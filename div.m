function out = div(vec,L)
    dx = compute_diff_ps(vec(:,:,1), 1, L(1));
    dy = compute_diff_ps(vec(:,:,2)', 1, L(2))';
    out = dx + dy;
end