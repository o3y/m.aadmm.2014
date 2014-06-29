function [p1new, p2new] = L2BallProjection(p1, p2)
    tmp = sqrt(p1.^2 + p2.^2);
    tmp(tmp < 1) = 1;
    p1new = p1 ./ tmp;
    p2new = p2 ./ tmp;
end
