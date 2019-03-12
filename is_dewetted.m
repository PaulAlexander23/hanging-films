function out = is_dewetted(H)
   
    out = any(any(H < 0));
    
end
