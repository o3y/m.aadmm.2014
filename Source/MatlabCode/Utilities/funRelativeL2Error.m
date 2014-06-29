function RelErr = funRelativeL2Error(x, x0)
    RelErr = norm(x(:) - x0(:)) / norm(x0(:));
end