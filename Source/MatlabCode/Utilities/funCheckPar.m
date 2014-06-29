function [par_val, par] =  funCheckPar(par, str_field, default_value)
% Check str_field in par. If field does not exist, assign default value to
% the field.
    if ~isfield(par, str_field)
        par.(str_field) = default_value;
    end
    par_val = par.(str_field);
end
