function [Etot, sEnergy] = funObjOverGraph(A, x, b, K)
% Calculates 1/2|Ax-b|^2 + ||Kx||_{4,1}
% where Kx = wReg * R(x)

    % Fidelity
    tmp = A * x - b;
    Efid = .5 * sum(tmp(:)' * tmp(:));
    
    % Regularity
    Kx = K * x;
    tmp = sum(sqrt(abs(Kx(1:4:end)).^2 + abs(Kx(2:4:end)).^2 ...
        + abs(Kx(3:4:end)).^2 + abs(Kx(4:4:end)).^2));
    Ereg = sum(tmp(:));
    
    % Total energy
    Etot = Efid + Ereg;
    
    sEnergy = [];
    sEnergy.Etot = Etot;
    sEnergy.Efid = Efid;
    sEnergy.Ereg = Ereg;
end
