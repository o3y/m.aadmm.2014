function h = funPlot(x, par)
    % Plot image x, and its residual (if available)
    h = check_par(par, 'hImage', []);
    xerr = check_par(par, 'Residual', []);
    iter = check_par(par, 'Iteration', []);
    sImageTitle = check_par(par, 'sImageTitle', 'Image');
    sResidualTitle = check_par(par, 'sResidualTitle', 'Residual');
    PauseTime = check_par(par, 'PauseTime', 0);
    cmap = check_par(par, 'ColorMap', 'gray');
    
    if ~isreal(x)
        x = real(x);
    end
    if ~isreal(xerr)
        xerr = abs(xerr);
    end
    
    if ~isempty(h) && ishandle(h)
        set(0,'CurrentFigure',h);
    else
        h = figure;
    end
    
    if ~isempty(iter)
        sImageTitle = ['Iteration ', num2str(iter), ', ', sImageTitle];
    end
    
    subplot(121);
    if isvector(x)
        plot(x, '.', 'MarkerSize', 20);
    else
        imagesc(x); 
        colormap(cmap); 
        axis image off;
    end
    title(sImageTitle);
    
    subplot(122);
    if isvector(xerr)
        plot(xerr);
    else
        imagesc(xerr); 
        colormap(gray); 
        axis image off;
    end
    title(sResidualTitle);            
    
    drawnow;
    pause(PauseTime);
end
