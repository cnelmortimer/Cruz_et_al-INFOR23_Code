function [ref, resName, resX, resVal] = mxParseConfs(inputStruct)
    fields = fieldnames(inputStruct);    
    resVal = -inf;
    focus = 1;
    
    for i=1:1:numel(fields)
        buffer = inputStruct.(fields{i});
        if buffer{3}(1) > resVal
           focus = i;
           resVal = buffer{3}(1);
        end
    end
    
    ref = fields{focus};
    chosen = inputStruct.(ref);    
    resName = chosen{1}{1};
    resX = chosen{2}(1,:);
    resVal = chosen{3}(1);
end
