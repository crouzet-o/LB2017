function [coords, theta] = lpoints(refpoint, distance, slope, dims)
    % Get the coordinates of one or several points on a line that is
    % defined by a reference point and a slope.
    % 
    % refpoint must have its coordinates in a 2D space
    % + for formants dimensions must be reversed!

    if isempty(dims)
        dims = [1, 2];
    end
    dim = size(refpoint);
    theta = atan(slope);
    
    X = refpoint(dims(1)) + (distance * cos(theta));
    Y = refpoint(dims(2)) + (distance * sin(theta));

    if (dim(2) > 1)
        coords = [X; Y];
    elseif (dim(1) > 1)
        coords = [X, Y];
    end
end