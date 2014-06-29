function funCropEdge(hFigure)
    ha = gca(hFigure);
    set(ha, 'Position', get(ha, 'OuterPosition') - ...
        get(ha, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
end