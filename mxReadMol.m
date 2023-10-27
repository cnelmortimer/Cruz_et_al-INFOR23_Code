function matrix = mxReadMol(fileName)
    file = fopen(fileName, 'r');
    if file<0
       error('File not found!'); 
    end
    
    fgetl(file); % Reading the '@<TRIPOS>MOLECULE'
    fgetl(file); % Reading the 'DB00...'
    numbers = fscanf(file, '%d', [1 5]);
    outNumber = numbers(1); % Accessing the first number of the array
    
    for i=1:1:5
        fgetl(file); % New line after fscanf, 'SMALL', 'USER_CHARGES', empty line, '@<TRIPOS>ATOM'
    end
    
    matrix = [];
    counter = 0;
    while ~feof(file)
        bufferLine = fgetl(file);
        if isequal(bufferLine, '@<TRIPOS>BOND') || counter==outNumber
            break;
        else
            pieces = split(bufferLine);
            pieces(1) = []; % Removing the previous break line
            matrix = [matrix; str2double(pieces{3}), str2double(pieces{4}), str2double(pieces{5})]; %#ok<AGROW>
            counter = counter + 1;
        end
    end
    fclose(file);
    if outNumber ~= size(matrix, 1)
       error('Inconsistent number of elements!'); 
    end
    matrix = matrix'; % Save every molecule in a column for faster access
end
