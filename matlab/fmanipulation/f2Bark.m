function output = f2Bark(input, inv)
    if nargin < 2
        inv = false;
    else
        inv = true;
    end
    if (~inv)
        output = (26.81 ./ (1 + 1960 ./ input )) - 0.53;
    else
        output = 1960 ./ (26.81 ./ (input + 0.53) - 1);
    end
end
