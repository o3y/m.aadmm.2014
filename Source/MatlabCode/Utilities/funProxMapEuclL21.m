function y = funProxMapEuclL21(y, dy)
    y = y + dy;
    tmp = sqrt(abs(y(:, :, 1)).^2 + abs(y(:, :, 2)).^2);
    tmp(tmp < 1) = 1;
    y(:, :, 1) = y(:, :, 1) ./ tmp;
    y(:, :, 2) = y(:, :, 2) ./ tmp;
end