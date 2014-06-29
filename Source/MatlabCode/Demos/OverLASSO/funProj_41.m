function w = funProj_41(w0, tau)
% Soft-thresholding for {k,1} norm
% This function solves the following optimization problem:
%   w = argmin F(w) + tau*|w - w0|^2/2
% where F(w) = |w|_{k,1}

    tmp = w0.*conj(w0);
    w0norm = sqrt(tmp(1:4:end) + tmp(2:4:end) + tmp(3:4:end) + tmp(4:4:end));
    weight = max(w0norm - 1/tau, 0) ./ max(w0norm, eps);
    w = w0;
    w(1:4:end) = w(1:4:end) .* weight;
    w(2:4:end) = w(2:4:end) .* weight;
    w(3:4:end) = w(3:4:end) .* weight;
    w(4:4:end) = w(4:4:end) .* weight;
end