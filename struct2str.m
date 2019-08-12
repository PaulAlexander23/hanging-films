function str = struct2str(struct)
    fields = fieldnames(struct);
    values = struct2cell(struct);
    str = "";
    for i = 1:length(fields)
        str = str + sprintf("-%s-%g", fields{i}, values{i}(end));
    end
end