function y = interp_over_nan(x)

    % Interpolates over NaNs and extrapolate on the sides by repeating the edge
    % values

    y = x;

    i = find(~isnan(y),1,'first');
    y(1:i-1) = y(i);

    i = find(~isnan(y),1,'last');
    y(i+1:end) = y(i);

    t = 1:length(x);
    s = ~isnan(y);

    y = interp1(t(s), y(s), t, 'linear');

end
