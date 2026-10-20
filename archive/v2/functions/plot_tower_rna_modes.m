function plot_tower_rna_modes(towerNodeData, ...
    modeShapeUndamped_global, naturalFreqUndamped_Hz, modeIndexToPlot)
%PLOT_TOWER_RNA_MODES 绘制选定模态的塔筒平动分量随高度的变化。
% 输入：完整振型包含塔底零自由度，每列一阶；阶次可为列向量。
% 输出：每阶一个图，分别展示 SS、FA 与竖直平动分量。
% 三个分量使用同一显示比例，以保留方向之间的相对大小。
% 转动分量仍保留在模型结果中；此图只展示平动，不用于模态定性分类。

height = towerNodeData.xyz(:, 3);
directionLabel = {'SS / uX', 'FA / uY', 'Vertical / uZ'};

for modeIdx = modeIndexToPlot(:).'
    nodeMode = reshape(modeShapeUndamped_global(:, modeIdx), 6, []).';
    translation = nodeMode(:, 1:3);
    displayScale = max(abs(translation(:)));
    if displayScale > 0
        translation = translation / displayScale;
    end

    figure('Name', sprintf('Tower + RNA mode %d', modeIdx));
    tiledlayout(1, 3);
    for directionIdx = 1:3
        nexttile;
        plot(translation(:, directionIdx), height, '-o', 'LineWidth', 1.2);
        xlabel(directionLabel{directionIdx});
        ylabel('Height (m)');
        xlim([-1.05, 1.05]);
        grid on;
    end
    sgtitle(sprintf('Mode %d: %.6f Hz (relative shape)', ...
        modeIdx, naturalFreqUndamped_Hz(modeIdx)));
end
end
