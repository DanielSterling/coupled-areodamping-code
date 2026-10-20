function modalOut = solve_tower_rna_undamped_modes( ...
    K_free, M_free, numModesRequested)
%SOLVE_TOWER_RNA_UNDAMPED_MODES 求约束后结构的无阻尼模态。
%
% 输入：K_free、M_free 为塔底约束后的刚度、质量矩阵；
%       numModesRequested 为请求的最低阶模态数。
% 输出：频率与角频率均为列向量；振型每列对应一阶模态。
%       振型按广义质量归一化，不表示真实位移幅值。
%       请求阶数超过模型自由度时，返回全部模态。
%
% 求解 K*phi = omega^2*M*phi，排序后令 phi'*M*phi=1。
% 当前只求无阻尼模态，不提前构造非比例阻尼或状态空间方程。

numFreeDof = size(K_free, 1);
numModesSolved = min(numModesRequested, numFreeDof);

if numModesSolved == numFreeDof
    [modeShapeUndamped, eigenvalueMatrix] = eig(full(K_free), full(M_free));
    eigSolveFlag = 0;
else
    eigOptions.issym = true;
    eigOptions.isreal = true;
    eigOptions.tol = 1e-10;
    eigOptions.maxit = 1000;
    eigOptions.p = min(numFreeDof, max(2*numModesSolved + 10, 30));
    [modeShapeUndamped, eigenvalueMatrix, eigSolveFlag] = ...
        eigs(K_free, M_free, numModesSolved, 0, eigOptions);
end

% 未收敛或非正特征值意味着本阶段结构结果不能用于后续计算。
if eigSolveFlag ~= 0
    error('无阻尼模态求解未收敛，eigs flag = %d。', eigSolveFlag);
end
eigenvalue = real(diag(eigenvalueMatrix));
if any(eigenvalue <= 0)
    error('固定塔底结构出现非正特征值，请检查结构矩阵和边界条件。');
end

[eigenvalue, modeSortIndex] = sort(eigenvalue, 'ascend');
modeShapeUndamped = real(modeShapeUndamped(:, modeSortIndex));

for modeIdx = 1:numModesSolved
    modeVector = modeShapeUndamped(:, modeIdx);
    modalMass = modeVector' * M_free * modeVector;
    modeShapeUndamped(:, modeIdx) = modeVector / sqrt(modalMass);
end

naturalOmegaUndamped_rad_s = sqrt(eigenvalue);
modalOut.numModesSolved = numModesSolved;
modalOut.eigSolveFlag = eigSolveFlag;
modalOut.modeShapeUndamped = modeShapeUndamped;
modalOut.naturalOmegaUndamped_rad_s = naturalOmegaUndamped_rad_s;
modalOut.naturalFreqUndamped_Hz = naturalOmegaUndamped_rad_s / (2*pi);
end
