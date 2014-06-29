function w = funThres21Norm(w0, tau)
% Soft-thresholding for {2,1} norm
% This function solves the following optimization problem:
%   w = argmin F(w) + tau*|w - w0|^2/2
% where F(w) = |w|_{2,1}

    w0norm = sqrt(sum(w0.*conj(w0), 3));
    tmp = max(w0norm - 1/tau, 0) ./ max(w0norm, eps);
    w = w0;
    w(:, :, 1) = w(:, :, 1) .* tmp;
    w(:, :, 2) = w(:, :, 2) .* tmp;
end