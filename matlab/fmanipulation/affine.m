function [slope, yint] = affine(A, B)

    slope = (B(2)-A(2)) / (B(1)-A(1));
    yint = B(2) - slope * B(1);

end
