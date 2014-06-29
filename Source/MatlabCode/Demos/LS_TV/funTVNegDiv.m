function Kty = funTVNegDiv(y, lambda, bPBC)
% The negative divergence operator: KTy = -div(y)
    [nRow, nCol, ~] = size(y);
    if bPBC
        % Periodic boundary condition
        Kty = diff([y(end,:,1); y(:,:,1)], 1, 1) + diff([y(:,end,2), y(:,:,2)], 1, 2);
        Kty = - Kty * lambda;
    else
        D1p = lambda * diff([zeros(1, nCol); y(1:end-1, :, 1); zeros(1, nCol)], 1, 1);
        D2p = lambda * diff([zeros(nRow, 1), y(:, 1:end-1, 2), zeros(nRow, 1)], 1, 2);
        Kty = - D1p - D2p;
    end        
end