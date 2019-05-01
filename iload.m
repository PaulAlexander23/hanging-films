function y0 = iload(~,filename)
    load(filename,'h');
    y0 = h;
end
