function [xag, etc] = funLS_TV_ADMM(A, b, wTV, par)
% L-ADMM and LP-ADMM

LipG = par.LipG;
LipK = par.LipK;
xsize = par.xsize;
xTrue = par.xTrue;
bVerbose = funCheckPar(par, 'bVerbose', true);
chi = double(funCheckPar(par, 'bPreconditioned', 0));
fhProjx = funCheckPar(par, 'fhProjx', @(x)(x));
MaxIter = funCheckPar(par, 'MaxIter', 100);
RhoEst = funCheckPar(par, 'RhoEstimate', 1);
bErgodic = funCheckPar(par, 'bErgodic', 0);

FKtKFt = wTV^2 * (abs(psf2otf([1,-1],xsize)).^2 + abs(psf2otf([1;-1],xsize)).^2);

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
w = zeros([xsize, 2]);
y = zeros([xsize, 2]);
Kx = zeros([xsize, 2]);
tStart = tic;

for t = 1:MaxIter
    % --------------------------------------
    % Main iteration
    % --------------------------------------
    
    % ------Middle step
    gradx = A'*(A*x - b);
    
    % ------x iteration (with or without preconditioning)
    if chi
        x = x - (funTVNegDiv(theta*(Kx - w) + y, wTV, 1) + gradx)/eta;
        x = fhProjx(x);
    else
        %   x_{t+1} = argmin theta_t*|Kx - Bw_t + b|^2/2 + <gradG, x> + <y_t, Kx> + eta_t*|x -x_t|^2/2
        %           = argmin theta_t*|Kx|^2/2 + <x, K'(theta_t*(b-Bw_t)+y_t) + gradG> + + eta_t*|x -x_t|^2/2
        RHS = funTVNegDiv(theta*w - y, wTV, 1) - gradx + eta*x;
        LHS = theta*FKtKFt + eta;
        x = ifft2(fft2(RHS) ./ LHS);
    end
    Kx = funTVGrad(x, wTV, 1);
    % ------w iteration
    w = funThres21Norm(Kx + y/tau, tau);
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