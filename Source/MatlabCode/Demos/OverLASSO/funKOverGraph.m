function [K, Kt] = funKOverGraph(nRow, nCol, wReg)
    if nargin < 3
        wReg = 1;
    end
    nnKRow = 1:(4*(nRow-1)*(nCol-1));
    nnKCol = zeros(4*(nRow-1)*(nCol-1), 1);
    ind = 1;
    for i = 1 : nCol-1
        for j = 1 : nRow-1
            nnKCol(ind:ind+3) = sub2ind([nRow, nCol], [j, j+1, j, j+1], [i, i, i+1, i+1]);
            ind = ind + 4;
        end
    end
    K = wReg * sparse(nnKRow, nnKCol, ones(1, (nRow-1)*(nCol-1)*4));
    Kt = K';
end