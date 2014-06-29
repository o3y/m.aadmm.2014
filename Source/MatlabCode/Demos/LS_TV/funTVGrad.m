function Kx = funTVGrad(x, lambda, bPBC)
% The finite difference operator: Kx = lambda * grad(x)
    [nRow, nCol] = size(x);
    Kx = zeros(nRow, nCol, 2);
    if bPBC
        % Periodic boundary condition
        Kx(:, :, 1) = diff([x; x(1, :)], 1, 1);
        Kx(:, :, 2) = diff([x, x(:, 1)], 1, 2);
        Kx = lambda * Kx;
    else
        Kx(:, :, 1) = lambda * [diff(x, 1, 1); zeros(1, nCol)];
        Kx(:, :, 2) = lambda * [diff(x, 1, 2), zeros(nRow, 1)];
    end
end

