% Script for comparison plots. Require variables xAAA, etcAAA, xBBB, etcBBB, ...
% The lAlg of variables should be saved as lAlg = {'AAA', 'BBB'}
if ~exist('scale', 'var')
    scale = [0, 1];
end
nResults = length(lAlg);

if bCompareObjVal
    hObjVal = figure;
end
if bCompareRelErr
    hRelErr = figure;
end
if bCompareImage
    hImg = figure;
    % ---------True Image------------
%     subplot(2, nResults+1, 1);
    subplot(1, 2, 1);
    imagesc(abs(xTrue), scale); colormap gray; axis equal off;
    title('True');
end
    
for i = 1:nResults
    sAlg = lAlg{i};
    eval(['tx = x', sAlg, ';']);
    eval(['tEtc = etc', sAlg, ';']);
    if ~isreal(tx)
        tx = real(tx);
    end
    tRes = abs(tx - xTrue);

    if bCompareImage
%         figure(hImg);
        % ---------Recovered Image-------
%         subplot(2, nResults + 1, i + 1);
        figure;
        subplot(1, 2, 1);
        imagesc(abs(tx), scale); colormap gray; axis equal off;
        title(lTitle{i}, 'Interpreter', 'none');
%         subplot(2, nResults + 1, nResults + i + 2);
        subplot(1, 2, 2);
        imagesc(tRes); colormap gray; axis equal off; % colorbar;
        title(sprintf('relerr=%g', funRelativeL2Error(tx, xTrue)));
    end
    if bCompareObjVal
        figure(hObjVal);

        if strcmp(sXLabel, 'Iteration')
            plot(tEtc.PrimalObjectiveValue, 'DisplayName', lTitle{i});        
            xlabel('Iteration');
        elseif strcmp(sXLabel, 'CPUTime')
            semilogy(tEtc.CPUTime, tEtc.PrimalObjectiveValue, 'DisplayName', lTitle{i});
            xlabel('CPU Time');
        end
        
        ylabel('Objective Value');
        hold all;
    end        
    if bCompareRelErr
        figure(hRelErr);
        
        if strcmp(sXLabel, 'Iteration')
            plot(tEtc.RelativeError, 'DisplayName', lTitle{i});
            xlabel('Iteration');
        elseif strcmp(sXLabel, 'CPUTime')
            plot(tEtc.CPUTime, tEtc.RelativeError, 'DisplayName', lTitle{i});
            xlabel('CPU Time');
        end
        
        ylabel('Relative Error to Ground Truth');
        hold all;
    end        
end

if bShowLS
    figure;
    subplot(2, 1, 1);
    imagesc(xLS); colormap gray; axis equal off; % colorbar;
    title(sprintf('relerr=%g', funRelativeL2Error(xLS, xTrue)));
    subplot(2, 1, 2);
    imagesc(abs(xLS - xTrue)); colormap gray; axis equal off; % colorbar;
end

if bCompareObjVal
    figure(hObjVal);
    if strcmp(sXLabel, 'Iteration')
        plot(TrueEnergy * ones(MaxIter, 1), 'DisplayName', 'True');
    elseif strcmp(sXLabel, 'CPUTime')
        plot([0, tEtc.CPUTime(end)], [TrueEnergy, TrueEnergy], 'DisplayName', 'True');
    end
    h = legend('show');
end
if bCompareRelErr
    figure(hRelErr);
    h = legend('show');
end
