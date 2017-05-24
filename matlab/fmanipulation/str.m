function s = str(a)

    if isnumeric(a)
        s = num2str(a);
    else
        s = a;
    end
end
