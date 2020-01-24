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
                funcStr = func2str(value);
                
                while contains(funcStr, '/')
                    funcStr = extractAfter(funcStr,'/');
                end
                
                str = str + sprintf("-%s-%s", fields{i}, funcStr);
            else
                str = str + sprintf("-%s-%g", fields{i}, value);
            end
        end
    end
end