function [xag, etc] = funLS_TV_AADMM_B(A, b, wTV, par)
% ALP-ADMM with bounded feasible sets and no backtracking

LipG = par.LipG;
LipK = par.LipK;
xsize = par.xsize;
xTrue = par.xTrue;
bVerbose = funCheckPar(par, 'bVerbose', true);
MaxIter = funCheckPar(par, 'MaxIter', 100);
RhoEst = funCheckPar(par, 'RhoEstimate', 1/LipK);
fhProjx = funCheckPar(par, 'fhProjx', @(x)(min(max(x, 0), 1)));

% --------------------------------------
% Initialization
% --------------------------------------
etc = [];
etc.CPUTime = nan(MaxIter, 1);
etc.RelativeError = nan(MaxIter, 1);
etc.PrimalObjectiveValue = nan(MaxIter, 1);

rho = RhoEst;
tau = RhoEst;
x = zeros(xsize);
xag = zeros(xsize);
w = zeros([xsize, 2]);
y = zeros([xsize, 2]);
Kx = zeros([xsize, 2]);
tStart = tic;

for t = 1:MaxIter
    % --------------------------------------
    % Main iteration
    % --------------------------------------
    % ------Stepsize parameters
    theta = (t-1)*RhoEst/t;
    alpha = 2/(t+1);
    eta = 2*LipG/t;
    
    % ------Middle step
    xmd = (1 - alpha) * xag + alpha * x;
    gradx = A'*(A*xmd - b);
    
    % ------x iteration
    x = x - (funTVNegDiv(theta*(Kx - w) + y, wTV, 1) + gradx)/eta;
    x = fhProjx(x);
    Kx = funTVGrad(x, wTV, 1);
    % ------w iteration
    w = funThres21Norm(Kx + y/tau, tau);
    % ------y iteration
    y = y - rho * (w - Kx);
    % ------Aggregate step
    xag = (1 - alpha) * xag + alpha * x;

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