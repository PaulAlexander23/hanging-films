function out = is_dewetted(H)
   
    out = all(all(H > 0))...
        * (abs(sum(sum(H - 1))) < 1e-5);
    
end