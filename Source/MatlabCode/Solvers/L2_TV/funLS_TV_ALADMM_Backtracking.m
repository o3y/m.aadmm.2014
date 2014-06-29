function [xag, etc] = funLS_TV_ALADMM_Backtracking(A, b, wTV, par)
% AL-ADMM with backtracking

xsize = par.xsize;
xTrue = par.xTrue;
bVerbose = funCheckPar(par, 'bVerbose', true);
MaxIter = funCheckPar(par, 'MaxIter', 100);
rhoEst = funCheckPar(par, 'RhoEstimate', 1);
L0 = funCheckPar(par, 'L0', 1);
Lmin = funCheckPar(par, 'Lmin', L0/100);
nu = funCheckPar(par, 'nu', (MaxIter-1)/Lmin);

FKtKFt = wTV^2 * (abs(psf2otf([1,-1],xsize)).^2 + abs(psf2otf([1;-1],xsize)).^2);

% --------------------------------------
% Initialization
% --------------------------------------
etc = [];
etc.CPUTime = nan(MaxIter, 1);
etc.RelativeError = nan(MaxIter, 1);
etc.PrimalObjectiveValue = nan(MaxIter, 1);
etc.Lt = nan(MaxIter, 1);
etc.LineSearchCount = nan(MaxIter, 1);
etc.alpha = nan(MaxIter, 1);
etc.Gamma = nan(MaxIter, 1);
etc.theta = nan(MaxIter, 1);
etc.tau = nan(MaxIter, 1);
etc.rho = nan(MaxIter, 1);

Gamma = L0;
L = L0;
x = zeros(xsize);
xag = zeros(xsize);
w = zeros([xsize, 2]);
y = zeros([xsize, 2]);
Resid = A*x - b;
gag = A' * Resid;
g = gag;

tStart = tic;
for t = 1:MaxIter
    % --------------------------------------
    % Main iteration
    % --------------------------------------
    % ------Variable updating
    Gamma_prev = Gamma;
    x_prev = x;
    Resid_prev = Resid;
    g_prev = g;
    gag_prev = gag;
    y_prev = y;
    w_prev = w;
    xag_prev = xag;
    
    bKeepSearch = 1;
    count = 1;
    while bKeepSearch
        alpha = (-Gamma_prev + sqrt(Gamma_prev^2 + 4*L*Gamma_prev)) / (2*L);
        if alpha>=1
            Gamma = L0;
        else
            Gamma = Gamma_prev * (1-alpha);
        end
        rho = alpha*rhoEst / (nu*Gamma);
        theta = nu*rhoEst*Gamma / alpha;
        tau = theta;
        eta = Gamma/alpha;

        % ------Middle step
        gmd = (1 - alpha) * gag_prev + alpha * g_prev;
        % ------x iteration
        RHS = funTVNegDiv(theta*w_prev - y_prev, wTV, 1) - gmd + eta*x_prev;
        LHS = theta*FKtKFt + eta;
        x = ifft2(fft2(RHS) ./ LHS);
        % Verify LipG
        Resid = A*x - b;
        Resid_diff = Resid(:) - Resid_prev(:);
        Resid_diff_norm = norm(Resid_diff);
        xdiff = x(:) - x_prev(:);
        xdiff_norm = norm(xdiff);
        if Resid_diff_norm > sqrt(L)*xdiff_norm
            count = count + 1;
            L = L*2;
        else
            etc.Lt(t) = L;
            etc.alpha(t) = alpha;
            etc.Gamma(t) = Gamma;
            etc.rho(t) = rho;
            etc.theta(t) = theta;
            etc.tau(t) = tau;
            etc.LineSearchCount(t) = count;
            etc.thetatEst(t) = nu;
            bKeepSearch = 0;
            g = A' * Resid;
            L = max(L/2, Lmin);
        end
    end
    % ------w iteration
    Kx = funTVGrad(x, wTV, 1);
    w = funThres21Norm(Kx + y_prev/tau, tau);
    % ------y iteration
    y = y_prev - rho * (w - Kx);
    % ------Aggregate step
    xag = (1 - alpha) * xag_prev + alpha * x;
    gag = (1 - alpha) * gag_prev + alpha * g;
    
    % --------------------------------------
    % Save CPU time
    % --------------------------------------
    etc.CPUTime(t) = toc(tStart);
    % --------------------------------------
    % Runtime outputs
    % --------------------------------------
    % Calculate the primal objective
    etc.PrimalObjectiveValue(t) = par.fhPrimalObjectiveValue(xag);
    % Calculate primal relative error to ground truth
    etc.RelativeError(t) = funRelativeL2Error(xag, xTrue);
    funPrintf(bVerbose, 't=%d,POBJ=%e,RelErr=%e\n', ...
        t, etc.PrimalObjectiveValue(t), etc.RelativeError(t));
end

end