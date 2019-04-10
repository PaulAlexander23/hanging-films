function y0 = iload(~,filename)
    load(filename,'h');
    y0 = reshape(h,64,64);
end