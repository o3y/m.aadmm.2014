% Demo for model sum_{i=1}^{nCH}|MFS_jx-f_j|^2/2 + wRegRescaled*TV(x), where M is
% the mask, F is the 2D Fourier transform matrix, S_j's are sensitivity
% maps, and f_j's are k-space datasets. The solution x of the model is the
% reconstructed image from partially parallel MRI.

%% Init
clear; clc; close all;

% General parameter
seed = 18;
MaxIter = 400;
RhoEstimate = 1;
sigma = 1e-4;
wReg = 1e-10;
beta = 10;
bSave = 0; % Flag for saving all results
sSaveDir = '../../../../Results/PPI';
sAlgRegexp = '(?!xTrue)^x.|^etc.|^par[^0].'; % Regular expression indicating algorithm results

% lAlg = {'BOSVS', 'ALADMML', 'ALPADMML'};
lAlg = {'ALADMML', 'ALPADMML'};
lMask = {'Linear', 'Poisson'};

% Parameters for output results
bCompareObjVal = 1;
bCompareImage = 0;
bCompareRelErr = 1;
bShowLS = 0;
% sXLabel = 'Iteration';
sXLabel = 'CPUTime';
lAlgAll = {'BOSVS', 'ALADMML', 'ALPADMML',};
lTitleAll = {'BOSVS', 'AL-ADMM w/ backtracking', 'ALP-ADMM w/ backtracking'};
[ism, ind] = ismember(lAlg, lAlgAll);
lAlg = lAlg(ism);
ind = ind(ism);
lTitle = lTitleAll(ind);

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
    diary([sRootDirName, '/logPPI.txt']);
    diary on;
    fprintf('%%-- %s --%% \r\n', datestr(now));
    fprintf('Results will be saved at %s\r\n', sRootDirName);
end
addpath(genpath('../../Solvers'));
addpath(genpath('../../Utilities'));
addpath(genpath('../../../../External/SBB_and_BOSVS'));
addpath('../../../../Data/PPI');
rng(seed);

%% Load problem
for iData = 1:2
    % for iData = 1
    for sMask = lMask
        load(sprintf('data%d.mat', iData));
        xTrue = u0;
        [nRow, nCol, nCh] = size(sense_map);
        switch sMask{:}
            case 'Poisson'
                %     Poisson mask
                mask = repmat(p, [1, 1, nCh]);
            case 'Linear'
                %     Linear mask
                mask = zeros(nRow, nCol, nCh);
                mask(1:6:end, :, :) = 1;
        end
        
        %     Non-orthonormal Fourier matrix; ||F||=sqrt(nRow*nCol), so we need to
        %     rescale wRegRescaled to wRegRescaled*nRow*nCol
        opA = @(x)(fft2(bsxfun(@times, x, sense_map)).*mask); % Matrix A
        opAt = @(y)(sum(ifft2(y.* mask) .* conj(sense_map), 3).*(nRow*nCol)); % Matrix A'
        LipG = max(max(abs(sum(sense_map, 3))))^2 * (nRow*nCol);
        wRegRescaled = wReg * nRow*nCol;
        sigmaRescaled = sigma * sqrt(nRow*nCol);
        
        LipK = sqrt(8) * wRegRescaled;
        
        objA = ClA_operator(opA, opAt);
        fTrue = opA(xTrue);
        Noise = mask .* (sigmaRescaled*(randn(nRow, nCol, nCh) + 1i*randn(nRow, nCol, nCh))/sqrt(2));
        fNoisy = fTrue + Noise;
        
        %% Necessary function handles, constants and parameters
        % Function handle for calculating energies
        fhPOBJ = @(x)(.5*sum(sum(sum(abs(objA*x-fNoisy).^2))) + sum(sum(sqrt(sum(abs(funTVGrad(x, wRegRescaled, 1)).^2,3)))));
        % Energy of ground truth
        TrueEnergy = fhPOBJ(xTrue);
        % Projection to the complex set {x: |x_i|<=1, i=1,...n}
        fhProjx = @(x)(x./max(abs(x),1));
        
        % Parameters for the uniform solver
        par = [];
        par.xsize = [nRow, nCol];
        par.Lmin = LipG / 10;
        par.L0 = nRow*nCol;
        par.M0 = LipK / 10;
        par.xTrue = xTrue;
        par.MaxIter = MaxIter;
        par.RhoEstimate = RhoEstimate;
        par.fhPrimalObjectiveValue = fhPOBJ;
        par0 = par;
        
        %     Save everything before running algorithms
        if bSave
            clear('-regexp', sAlgRegexp);
            sFilename = [sRootDirName, '/PPI', num2str(iData), sMask{:}, '.mat'];
            save(sFilename);
        end
        
        %% ALP-ADMM with backtracking
        sDescription = 'ALPADMML';
        sFunctionName = 'funLS_TV_ALPADMM_Backtracking';
        if any(strcmp(sDescription, lAlg))
            par = par0;
            par.fhProjx = fhProjx;
            
            eval(['par', sDescription, ' = par;']);
            fprintf('TVL2: %s algorithm is running...\n', sDescription);
            tic;
            [x, etc] = feval(sFunctionName, objA, fNoisy, wRegRescaled, par);
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
            [x, etc] = feval(sFunctionName, objA, fNoisy, wRegRescaled, par);
            toc;
            eval(['x', sDescription, '=x;']);
            eval(['etc', sDescription, '=etc;']);
        end
        

        %% Show results
        if bSave
            fprintf('Saving algorithm results...\r\n');
            save(sFilename, '-regexp', sAlgRegexp, '-append');
        end
        % Plot
        scriptComparisonPlot;
        % Save
        if bCompareObjVal
            figure(hObjVal);
            switch sMask{:}
                case 'Linear'
%                     ylim([TrueEnergy/2, (TrueEnergy+eps)*400]);
                case 'Poisson'
%                     ylim([TrueEnergy/2, (TrueEnergy+eps)*5]);
            end
            if bSave
                funCropEdge(hObjVal);
                print(hObjVal, '-depsc', ...
                    sprintf('%s/PPI%d%s_ObjVal_v_%s', sRootDirName, iData, sMask{:}, sXLabel));
            end
        end
        if bCompareRelErr
            figure(hRelErr);
            switch sMask{:}
                case 'Linear'
                    ylim([0, .6]);
                case 'Poisson'
                    ylim([0, .4]);
            end
            if bSave
                funCropEdge(hRelErr);
                print(hRelErr, '-depsc', ...
                    sprintf('%s/PPI%d%s_RelErr_v_%s', sRootDirName, iData, sMask{:}, sXLabel));
            end
        end
    end
end

if bSave
    close all;
    diary off;
    set(0, 'DefaultFigureVisible', 'on')
end

