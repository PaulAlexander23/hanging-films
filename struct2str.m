function str = struct2str(struct)
    fields = fieldnames(struct);
    values = struct2cell(struct);
    
    str = "";
    for i = 1:length(fields)
        value = values{i};
        if ~isempty(value)
            
            if isa(value, 'struct')
                str = str + struct2str(value);
            elseif isa(value, 'string') || isa(value, 'char')
                str = str + sprintf("-%s-%s", fields{i}, value);
            elseif isa(value, 'function_handle')
                 str = str + sprintf("-%s-%s", fields{i}, func2str(value));
            else
                str = str + sprintf("-%s-%g", fields{i}, value);
            end
        end
    end
end