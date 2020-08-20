function y0 = iload(~, filename)
    
    load(filename,'h');
    if ~exist('h', 'var') 
        load(filename,'solution');

        h = solution.y(:,:,end);
    end

    y0 = h;
end
