function [xag, etc] = funLS_TV_ALPADMM_Backtracking(A, b, wTV, par)
% ALP-ADMM with backtracking

xsize = par.xsize;
xTrue = par.xTrue;
bVerbose = funCheckPar(par, 'bVerbose', true);
MaxIter = funCheckPar(par, 'MaxIter', 100);
rhoEst = funCheckPar(par, 'RhoEstimate', 1);
L0 = funCheckPar(par, 'L0', 1);
Lmin = funCheckPar(par, 'Lmin', L0/100);
M0 = funCheckPar(par, 'M0', 1e-8);
fhProjx = funCheckPar(par, 'fhProjx', @(x)(min(max(x, 0), 1)));

% --------------------------------------
% Initialization
% --------------------------------------
etc = [];
etc.CPUTime = nan(MaxIter, 1);
etc.RelativeError = nan(MaxIter, 1);
etc.PrimalObjectiveValue = nan(MaxIter, 1);
etc.Lt = nan(MaxIter, 1);
etc.Mt = nan(MaxIter, 1);
etc.LineSearchCount = nan(MaxIter, 1);
etc.alpha = nan(MaxIter, 1);
etc.Gamma = nan(MaxIter, 1);
etc.theta = nan(MaxIter, 1);
etc.tau = nan(MaxIter, 1);
etc.rho = nan(MaxIter, 1);
etc.nu = nan(MaxIter, 1);

Gamma = L0;
L = L0;
M = M0;
nu = -Inf;
x = zeros(xsize);
xag = zeros(xsize);
Kx = zeros([xsize, 2]);
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
    nu_prev = nu;
    x_prev = x;
    Resid_prev = Resid;
    g_prev = g;
    gag_prev = gag;
    y_prev = y;
    w_prev = w;
    xag_prev = xag;
    Kx_prev = Kx;
    
    bKeepSearch = 1;
    count = 0;
    while bKeepSearch
        bKeepSearch = 0;
        count = count + 1;
        alpha = (-Gamma_prev + sqrt(Gamma_prev^2 + 4*L*Gamma_prev)) / (2*L);
        if alpha>=1
            Gamma = L0;
        else
            Gamma = Gamma_prev * (1-alpha);
        end
        nu = max(nu_prev, alpha/Gamma);
        rho = alpha*rhoEst / (nu*Gamma);
        theta = nu*rhoEst*Gamma / alpha;
        tau = theta;
        eta = Gamma/alpha + theta*M^2;

        % ------Middle step
        gmd = (1 - alpha) * gag_prev + alpha * g_prev;
        % ------x iteration
        x = x_prev - (funTVNegDiv(theta*(Kx_prev - w_prev) + y_prev, wTV, 1) + gmd)/eta;
        x = fhProjx(x);
        Kx = funTVGrad(x, wTV, 1);
        % Verify LipG and LipK
        Resid = A*x - b;
        Resid_diff = Resid(:) - Resid_prev(:);
        Resid_diff_norm = norm(Resid_diff);
        xdiff = x(:) - x_prev(:);
        xdiff_norm = norm(xdiff);
        if Resid_diff_norm > sqrt(L)*xdiff_norm
            L = L*2;
            bKeepSearch = 1;
        end
        if norm(Kx(:) - Kx_prev(:)) > M*xdiff_norm
            M = M*2;
            bKeepSearch = 1;
        end
    end
    etc.Lt(t) = L;
    etc.Mt(t) = M;
    etc.alpha(t) = alpha;
    etc.Gamma(t) = Gamma;
    etc.rho(t) = rho;
    etc.theta(t) = theta;
    etc.tau(t) = tau;
    etc.nu(t) = nu;
    etc.LineSearchCount(t) = count;
    g = A' * Resid;
    L = max(L/2, Lmin);
    % ------w iteration
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