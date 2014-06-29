% Demo for compressive sensing total-variation image reconstruction
%   min_x .5*|Ax-b|^2 + lambda*||Dx||_{2,1}
% where ||Dx||_{2,1} is discrete total-variation.

%% Init
clear; clc; close all;

% General parameter
seed = 18;
SampleRate = .3;
sigma = 1e-3;
MaxIter = 300;
wReg = 1e-3;
lProblem = {'Gaussian', 'Bernoulli'};
bBoxConstraint = 1;
xLB = 0;
xUB = 1;
lAlg = {'LADMM_U', 'LPADMM_Box', 'ALADMM_U', 'ALPADMM_Box', 'ALPADMM_U'};

bSave = 0; % Flag for saving all results
sSaveDir = '../../../../Results/CS';

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
    diary([sRootDirName, '/logCS.txt']);
    diary on;
    fprintf('%%-- %s --%% \r\n', datestr(now));
    fprintf('Results will be saved at %s\r\n', sRootDirName);
end
addpath(genpath('../../Solvers'));
addpath(genpath('../../Utilities'));
rng(seed);

for sProblem = lProblem
    %% Generate problem
    fprintf('Generating %s instance...\n', sProblem{:});
    xTrue = phantom(64);
    xTrue = min(max(xTrue, xLB), xUB);
    [nRow, nCol] = size(xTrue);
    nSample = ceil(SampleRate * nRow * nCol);
    switch sProblem{:}
        case 'Gaussian'
            A = randn(nSample, nRow*nCol) / sqrt(nSample);
        case 'Bernoulli'
            A = (1-(rand(nSample, nRow*nCol)<.5)*2) / sqrt(nSample);
        otherwise
            error('Unknown problem.')
    end
    fTrue = A * xTrue(:);
    Noise = randn(nSample, 1) * sigma;
    fNoisy = fTrue + Noise;
    objA = ClA_operator(@(x)(A*x(:)), @(f)(reshape((f'*A)', [nRow, nCol])));
    LipG = max(eig(A' * A));
    LipK = sqrt(8) * wReg;
    fprintf('L_G=%g, ||K||=%g\r\n', LipG, LipK);
    
    %% Necessary function handles, constants and parameters
    % Function handle for primal objective value
    fhPOBJ = @(x)(.5*sum(sum(sum(abs(objA*x-fNoisy).^2))) + sum(sum(sqrt(sum(abs(funTVGrad(x, wReg, 1)).^2,3)))));
    % Primal objective value of ground truth
    TrueEnergy = fhPOBJ(xTrue);
    % Projection to the box [xLB, xUB]
    fhProjx = @(x)(min(max(x, xLB), xUB));
    
    % Parameters for the uniform solver
    par = [];
    par.xsize = [nRow, nCol];
    par.LipG = LipG;
    par.LipK = LipK;
    par.Lmin = LipG / 10;
    par.L0 = par.Lmin;
    par.M0 = LipK / 10;
    par.xTrue = xTrue;
    par.MaxIter = MaxIter;
    par.fhPrimalObjectiveValue = fhPOBJ;
    par.RhoEstimate = 1/LipK;
    par0 = par;
    
    %% ALP-ADMM with backtracking
    sDescription = 'ALPADMML';
    sFunctionName = 'funLS_TV_ALPADMM_Backtracking';
    if any(strcmp(sDescription, lAlg))
        par = par0;
        par.fhProjx = fhProjx;
        
        eval(['par', sDescription, ' = par;']);
        fprintf('TVL2: %s algorithm is running...\n', sDescription);
        tic;
        [x, etc] = feval(sFunctionName, objA, fNoisy, wReg, par);
        toc;
        eval(['x', sDescription, '=x;']);
        eval(['etc', sDescription, '=etc;']);
    end
    
    %% AL-ADMM with backtracking
    sDescription = 'ALADMML';
    sFunctionName = 'funLS_TV_ALADMM_Backtracking';
    if any(strcmp(sDescription, lAlg))
        par = par0;
        
        eval(['par', sDescription, ' = par;']);
        fprintf('TVL2: %s algorithm is running...\n', sDescription);
        tic;
        [x, etc] = feval(sFunctionName, objA, fNoisy, wReg, par);
        toc;
        eval(['x', sDescription, '=x;']);
        eval(['etc', sDescription, '=etc;']);
    end
    %% ALP-ADMM, X=Box(xLB, xUB)
    sDescription = 'ALPADMM_Box';
    sFunctionName = 'funLS_TV_AADMM_B';
    if any(strcmp(sDescription, lAlg))
        par = par0;
        par.fhProjx = fhProjx;
        
        eval(['par', sDescription, ' = par;']);
        fprintf('TVL2: %s algorithm is running...\n', sDescription);
        tic;
        [x, etc] = feval(sFunctionName, objA, fNoisy, wReg, par);
        toc;
        eval(['x', sDescription, '=x;']);
        eval(['etc', sDescription, '=etc;']);
    end
    
    %% L-ADMM, X=R^n
    sDescription = 'LADMM_U';
    sFunctionName = 'funLS_TV_ADMM';
    if any(strcmp(sDescription, lAlg))
        par = par0;
        par.bPreconditioned = 0;
        
        eval(['par', sDescription, ' = par;']);
        fprintf('TVL2: %s algorithm is running...\n', sDescription);
        tic;
        [x, etc] = feval(sFunctionName, objA, fNoisy, wReg, par);
        toc;
        eval(['x', sDescription, '=x;']);
        eval(['etc', sDescription, '=etc;']);
    end
    
    %% LP-ADMM, X=Box(xLB, xUB)
    sDescription = 'LPADMM_Box';
    sFunctionName = 'funLS_TV_ADMM';
    if any(strcmp(sDescription, lAlg))
        par = par0;
        par.bPreconditioned = 1;
        par.fhProjx = fhProjx;
        
        eval(['par', sDescription, ' = par;']);
        fprintf('TVL2: %s algorithm is running...\n', sDescription);
        tic;
        [x, etc] = feval(sFunctionName, objA, fNoisy, wReg, par);
        toc;
        eval(['x', sDescription, '=x;']);
        eval(['etc', sDescription, '=etc;']);
    end
    
    %% AL-ADMM, X=R^n
    sDescription = 'ALADMM_U';
    sFunctionName = 'funLS_TV_AADMM_U';
    if any(strcmp(sDescription, lAlg))
        par = par0;
        par.bPreconditioned = 0;
        
        eval(['par', sDescription, ' = par;']);
        fprintf('TVL2: %s algorithm is running...\n', sDescription);
        tic;
        [x, etc] = feval(sFunctionName, objA, fNoisy, wReg, par);
        toc;
        eval(['x', sDescription, '=x;']);
        eval(['etc', sDescription, '=etc;']);
    end
    
    %% ALP-ADMM, X=R^n
    sDescription = 'ALPADMM_U';
    sFunctionName = 'funLS_TV_AADMM_U';
    if any(strcmp(sDescription, lAlg))
        par = par0;
        par.bPreconditioned = 1;
        
        eval(['par', sDescription, ' = par;']);
        fprintf('TVL2: %s algorithm is running...\n', sDescription);
        tic;
        [x, etc] = feval(sFunctionName, objA, fNoisy, wReg, par);
        toc;
        eval(['x', sDescription, '=x;']);
        eval(['etc', sDescription, '=etc;']);
    end
    
    %% Show results
    % Parameters for comparing algorithms
    bCompareObjVal = 1;
    bCompareImage = 0;
    bCompareRelErr = 1;
    bShowLS = 0;
    % sXLabel = 'Iteration';
    sXLabel = 'CPUTime';
    lAlgAll = {'LADMM_U', 'LPADMM_Box', 'ALADMM_U', 'ALPADMM_U', 'ALPADMM_Box', ...
        'ALADMML', 'ALPADMML'};
    lTitleAll = {'L-ADMM, X=R^n', 'LP-ADMM, bounded X',...
        'AL-ADMM, X=R^n', 'ALP-ADMM, X=R^n', 'ALP-ADMM, bounded X',...
        'AL-ADMM w/ backtracking', 'ALP-ADMM w/ backtracking'};
    
    [ism, ind] = ismember(lAlg, lAlgAll);
    lAlg = lAlg(ism);
    ind = ind(ism);
    lTitle = lTitleAll(ind);
    
    scriptComparisonPlot;
    
    if bSave
        fprintf('Saving figures...\r\n');
        if bCompareObjVal
            figure(hObjVal);
            ylim([TrueEnergy/2, (TrueEnergy+eps)*3]);
            funCropEdge(hObjVal);
            print(hObjVal, '-depsc',...
                sprintf('%s/%s_ObjVal_v_%s.eps', sRootDirName, sProblem{:}, sXLabel));
        end
        if bCompareRelErr
            figure(hRelErr);
            ylim([0,1]);
            funCropEdge(hRelErr);
            print(hRelErr, '-depsc',...
                sprintf('%s/%s_RelErr_v_%s.eps', sRootDirName, sProblem{:}, sXLabel));
        end
    end
end

if bSave
    close all;
    set(0, 'DefaultFigureVisible', 'on')
    diary off;
end