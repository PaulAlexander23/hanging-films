function E = energy(H,L)

    N = size(H);

    S = L'./N(1:2);

    E = sum(sum((H - 1).^2*prod(S)));
    
    E = squeeze(E);
end
