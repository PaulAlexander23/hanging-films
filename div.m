function out = div(vec,L)
    vec(:,:,1) = diff_ps(vec(:,:,1), 1, L(1));
    vec(:,:,2) = diff_ps(vec(:,:,2)', 1, L(2))';
    out = vec(:,:,1) + vec(:,:,2);
end