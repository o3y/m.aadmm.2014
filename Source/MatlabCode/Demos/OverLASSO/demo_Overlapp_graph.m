% function demo_Overlapp_graph
% Demo for model
%   min_x 1/2|Ax-b|^2 + lambda*R(x),
% where R is a regularization function, and A is randomly generated.

% The dataset in this demo is generated according to section 9.3 of the
% following paper:
% Jacob, Laurent, Guillaume Obozinski, and Jean-Philippe Vert. "Group lasso
% with overlap and graph lasso." In Proceedings of the 26th Annual
% International Conference on Machine Learning, pp. 433-440. ACM, 2009.

%% Init
clear;
clc; close all;

%% Parameters
% General parameter
seed = 5; nRow = 64; nCol = 64;
sigma = .01;
nSample = floor(.5*nRow*nCol);
wReg = 1;
MaxIter = 300;

bSave = 0; % Flag for saving all results
sSaveDir = '../../../../Results/OverLASSO';
sProblem = '4cycle';
lAlg = {'LADMM', 'LPADMM', 'ALADMM', 'ALPADMM'};


%% Generate problem
if bSave
    set(0, 'DefaultFigureVisible', 'off');
    bLocked = 1;
    while bLocked
        sRootDirName = [sSaveDir, '/', num2str(now, '%.4f')];
        if ~exist(sRootDirName, 'dir')
            mkdir(sRootDirName);
            bLocked = 0;
        end
    end
    diary([sRootDirName, '/logOverLASSO.txt']);
    diary on;
    fprintf('%%-- %s --%% \r\n', datestr(now));
    fprintf('Results will be saved at %s\r\n', sRootDirName);
end
addpath(genpath('../../Solvers'));
addpath(genpath('../../Utilities'));
addpath(genpath('../../../../Data'));
rng(seed, 'twister');

switch sProblem
    case '4cycle'
        % Generate domain of x
        xSupp = false(nRow, nCol);
        for i = 1:3
            ulRow = randi(nRow);
            ulCol = randi(nCol);
            lrRow = ulRow + randi(nRow - ulRow) - 1;
            lrCol = ulCol + randi(nCol - ulCol) - 1;
            xSupp(ulRow:lrRow, ulCol:lrCol) = 1;
        end
        % Generate A
        A = randn(nSample, nRow*nCol);
        
        % Generate ground truth: random values in five randomly located
        % rectangles
        xTrue = zeros(nRow, nCol);
        nSupp = nnz(xSupp);
        xTrue(xSupp) = randn(nSupp, 1);
        % Generate measurements
        fTrue = A * xTrue(:);
        Noise = randn(nSample, 1) * sigma;
        fNoisy = fTrue + Noise;
        LipG = max(eig(A'*A));
end
% Operators K, Kt
K = funKOverGraph(nRow, nCol, wReg);
KtK = spdiags(K' * K);
LipK = sqrt(max(abs(KtK))); % Should be 2*wReg

%% Generate necessary operators & constants
% Function handle for calculating energies
fhPOBJ = @(x)(funObjOverGraph(A, x, fNoisy, K));
% Calculate objective value of ground truth
TrueEnergy = fhPOBJ(xTrue(:));

par = [];
par.xsize = [nRow*nCol, 1];
par.wsize = [4*(nRow-1)*(nCol-1), 1];
par.LipG = LipG;
par.LipK = LipK;
par.xTrue = xTrue(:);
par.MaxIter = MaxIter;
par.fhPrimalObjectiveValue = fhPOBJ;
par.RhoEstimate = 1/LipK;
par0 = par;

%% L-ADMM
sDescription = 'LADMM';
sFunctionName = 'funLS_DiagKtK_ADMM';
if any(strcmp(sDescription, lAlg))
    par = par0;
    par.bPreconditioned = 0;
    
    eval(['par', sDescription, ' = par;']);
    fprintf('TVL2: %s algorithm is running...\n', sDescription);
    tic;
    [x, etc] = feval(sFunctionName, A, fNoisy, K, KtK, par);
    x = reshape(x, [nRow, nCol]);
    toc;
    eval(['x', sDescription, '=x;']);
    eval(['etc', sDescription, '=etc;']);
end

%% LP-ADMM
sDescription = 'LPADMM';
sFunctionName = 'funLS_DiagKtK_ADMM';
if any(strcmp(sDescription, lAlg))
    par = par0;
    par.bPreconditioned = 1;
    
    eval(['par', sDescription, ' = par;']);
    fprintf('TVL2: %s algorithm is running...\n', sDescription);
    tic;
    [x, etc] = feval(sFunctionName, A, fNoisy, K, KtK, par);
    x = reshape(x, [nRow, nCol]);
    toc;
    eval(['x', sDescription, '=x;']);
    eval(['etc', sDescription, '=etc;']);
end

%% AL-ADMM
sDescription = 'ALADMM';
sFunctionName = 'funLS_DiagKtK_AADMM_U';
if any(strcmp(sDescription, lAlg))
    par = par0;
    par.bPreconditioned = 0;
    
    eval(['par', sDescription, ' = par;']);
    fprintf('TVL2: %s algorithm is running...\n', sDescription);
    tic;
    [x, etc] = feval(sFunctionName, A, fNoisy, K, KtK, par);
    x = reshape(x, [nRow, nCol]);
    toc;
    eval(['x', sDescription, '=x;']);
    eval(['etc', sDescription, '=etc;']);
end

%% ALP-ADMM
sDescription = 'ALPADMM';
sFunctionName = 'funLS_DiagKtK_AADMM_U';
if any(strcmp(sDescription, lAlg))
    par = par0;
    par.bPreconditioned = 1;
    
    eval(['par', sDescription, ' = par;']);
    fprintf('TVL2: %s algorithm is running...\n', sDescription);
    tic;
    [x, etc] = feval(sFunctionName, A, fNoisy, K, KtK, par);
    x = reshape(x, [nRow, nCol]);
    toc;
    eval(['x', sDescription, '=x;']);
    eval(['etc', sDescription, '=etc;']);
end

%% Show result;
% Parameters for comparing algorithms
bCompareObjVal = 1;
bCompareImage = 0;
bCompareRelErr = 1;
bShowLS = 0;
% sXLabel = 'Iteration';
sXLabel = 'CPUTime';

lAlgAll = {'LADMM', 'LPADMM', 'ALADMM', 'ALPADMM'};
lTitleAll = {'L-ADMM', 'LP-ADMM', 'AL-ADMM', 'ALP-ADMM'};

[ism, ind] = ismember(lAlg, lAlgAll);
lAlg = lAlg(ism);
ind = ind(ism);
lTitle = lTitleAll(ind);

% Plot
scriptComparisonPlot;
% Save
if bSave
    fprintf('Saving figures...\r\n');
    if bCompareObjVal
        figure(hObjVal);
        ylim([TrueEnergy/2, (TrueEnergy+eps)*3]);
        funCropEdge(hObjVal);
        print(hObjVal, '-depsc', ...
            sprintf('%s/OverLASSO_%s_ObjVal_v_%s.eps', sRootDirName, sProblem, sXLabel));
    end
    if bCompareRelErr
        figure(hRelErr);
        ylim([0,1]);
        funCropEdge(hRelErr);
        print(hRelErr, '-depsc', ...
            sprintf('%s/OverLASSO_%s_RelErr_v_%s.eps', sRootDirName, sProblem, sXLabel));
    end
    close all;
    diary off;
    set(0, 'DefaultFigureVisible', 'on')
end
