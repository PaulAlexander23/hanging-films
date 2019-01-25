function V = volume(H,L)

    N = size(H);

    S = L./N(1:2);

    V = sum(sum(H*prod(S)));
    
end
