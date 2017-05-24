function s = make_fname(filename, p)

    s = {};
    p = orderfields(p);
    fields = fieldnames(p);
    for f=fields
        s{end+1} = sprintf('%s=%s', f{1}, str(p.(f{1})));
    end
    
    [d, base, ~] = fileparts(filename);

    %s = fullfile(d, [base, '_', implode(s, '_')]);
    s = fullfile(d, [base, '_', [sprintf('%s_', s{1:end-1}), s{end}]]);
end

%======================================================================
