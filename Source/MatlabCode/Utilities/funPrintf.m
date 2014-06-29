function funPrintf(bVerbose, varargin)
    if bVerbose
        fprintf(varargin{:});
    end
end