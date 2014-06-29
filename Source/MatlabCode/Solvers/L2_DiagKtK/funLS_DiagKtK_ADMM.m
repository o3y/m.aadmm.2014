function [xag, etc] = funLS_DiagKtK_ADMM(A, b, K, KtK, par)

LipG = par.LipG;
LipK = par.LipK;
xsize = par.xsize;
wsize = par.wsize;
xTrue = par.xTrue;
bVerbose = funCheckPar(par, 'bVerbose', true);
chi = double(funCheckPar(par, 'bPreconditioned', 0));
MaxIter = funCheckPar(par, 'MaxIter', 100);
RhoEst = funCheckPar(par, 'RhoEstimate', 1);
bErgodic = funCheckPar(par, 'bErgodic', 0);
Kt = K';

% --------------------------------------
% Initialization
% --------------------------------------
etc = [];
etc.CPUTime = nan(MaxIter, 1);
etc.RelativeError = nan(MaxIter, 1);
etc.PrimalObjectiveValue = nan(MaxIter, 1);

eta = LipG + chi * RhoEst * LipK^2;
rho = RhoEst;
tau = rho;
theta = rho;

x = zeros(xsize);
xag = zeros(xsize);
w = zeros(wsize);
y = zeros(wsize);
Kx = zeros(wsize);
tStart = tic;

for t = 1:MaxIter
    % --------------------------------------
    % Main iteration
    % --------------------------------------
    
    % ------Middle step
    gradx = A'*(A*x - b);
    
    % ------x iteration (with or without preconditioning)
    if chi
        x = x - (Kt*(theta*(Kx - w) + y) + gradx)/eta;
    else
        %   x_{t+1} = argmin theta_t*|Kx - Bw_t + b|^2/2 + <gradG, x> + <y_t, Kx> + eta_t*|x -x_t|^2/2
        %           = argmin theta_t*|Kx|^2/2 + <x, K'(theta_t*(b-Bw_t)+y_t) + gradG> + + eta_t*|x -x_t|^2/2
        x = (Kt*(theta*w - y) - gradx + eta*x) ./ (theta*KtK + eta);
    end
    Kx = K*x;
    % ------w iteration
    w = funProj_41(Kx + y/tau, tau);
    % ------y iteration
    y = y - rho * (w - Kx);
    % ------Aggregate step
    if bErgodic
        xag = (xag * (t-1) + xnew) / t;
    else
        xag = x;
    end

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