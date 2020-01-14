function str = struct2str(struct)
    fields = fieldnames(struct);
    values = struct2cell(struct);
    
    str = "";
    for i = 1:length(fields)
        value = values{i};
        if ~isempty(value)
            
            if isa(value, 'function_handle')
                value = func2str(value);
            end
            if isa(value, 'struct')
                str = str + struct2str(value);
            else
                str = str + sprintf("-%s-%s", fields{i}, num2str(value));
            end
        end
    end
end