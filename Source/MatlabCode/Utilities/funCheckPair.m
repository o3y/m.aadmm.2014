% FUNCHECKPAIR   Check and assign default value for parameter pairs
%   [bFun, fhFun] = funCheckPair(SPar, s_bFun, s_fhFun) checks whether both
%   boolean flag bFun and function handle fhFun exist in paramter structure
%   SPar. If any of them does not exist in SPar, bFun is assigned the
%   default value (false).

function [bFun, fhFun] = funCheckPair(SPar, s_bFun, s_fhFun)
    fhFun = funCheckPar(SPar, s_fhFun, []);
    bFun = funCheckPar(SPar, s_bFun, false) && ~isempty(fhFun);
end