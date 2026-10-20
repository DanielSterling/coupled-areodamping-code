function towerModel = assemble_tower_rna_structure( ...
    towerNodeData, towerElemData, towerConn, massMatrixRNA_atTowerTop)
%ASSEMBLE_TOWER_RNA_STRUCTURE 装配柔性塔筒与塔顶 RNA 集中质量。
%
% 输入：
%   towerNodeData.xyz       nNode x 3 塔筒节点全局坐标。
%   towerElemData          原塔筒预处理函数生成的梁单元属性。
%   towerConn              nElem x 2 单元连接关系。
%   massMatrixRNA_atTowerTop 6 x 6，RNA 关于塔顶的全局质量矩阵。
%     自由度顺序必须为 [uX uY uZ thetaX thetaY thetaZ]。
%
% 输出 towerModel：
%   K_global、M_global     6*nNode 阶总体结构矩阵。
%   M_tower_global         单独的塔筒质量矩阵。
%   M_RNA_global           RNA 在总体自由度中的质量矩阵。
%   towerTopSelectionMatrix 6 x (6*nNode)，塔顶运动提取矩阵。
%   其余字段              节点、单元、自由度映射信息。
%
% 物理步骤：塔筒单元积分与装配 -> 塔顶刚体动能等效。
% 对于 RNA，T = 0.5*qdot_top'*M6*qdot_top，q_top=P*q，
% 所以 RNA 对总体质量矩阵的贡献为 P'*M6*P。
% M6 已包含参考点偏置的耦合项，不能再次平移或拆分重复装配。
% RNA 的变形忽略，不需要额外轮毂节点和高刚度刚臂单元。

dofPerNode = 6;
numTowerNode = size(towerNodeData.xyz, 1);
numTowerElem = numel(towerElemData);
numDofTotal = dofPerNode * numTowerNode;

nodeBaseDof1 = (towerConn(:, 1) - 1) * dofPerNode;
nodeBaseDof2 = (towerConn(:, 2) - 1) * dofPerNode;
elemDofMap = [nodeBaseDof1 + (1:dofPerNode), ...
              nodeBaseDof2 + (1:dofPerNode)];

K_global = sparse(numDofTotal, numDofTotal);
M_tower_global = sparse(numDofTotal, numDofTotal);

for towerElemIdx = 1:numTowerElem
    [K_elem_global, M_elem_global] = ...
        beam3d_tower_TIM_station(towerElemData(towerElemIdx));
    elemDof = elemDofMap(towerElemIdx, :);

    K_global(elemDof, elemDof) = ...
        K_global(elemDof, elemDof) + K_elem_global;
    M_tower_global(elemDof, elemDof) = ...
        M_tower_global(elemDof, elemDof) + M_elem_global;
end

towerTopNodeId = numTowerNode;
towerTopNodeDof = ((towerTopNodeId - 1)*dofPerNode + (1:dofPerNode)).';
towerTopSelectionMatrix = sparse( ...
    (1:dofPerNode).', towerTopNodeDof, ones(dofPerNode, 1), ...
    dofPerNode, numDofTotal);

M_RNA_global = towerTopSelectionMatrix' * ...
    sparse(massMatrixRNA_atTowerTop) * towerTopSelectionMatrix;
M_global = M_tower_global + M_RNA_global;

towerModel = struct();
towerModel.numTowerElem = numTowerElem;
towerModel.numTowerNode = numTowerNode;
towerModel.numDofTotal = numDofTotal;
towerModel.dofPerNode = dofPerNode;
towerModel.globalConn = towerConn;
towerModel.elemDofMap = elemDofMap;
towerModel.towerTopNodeId = towerTopNodeId;
towerModel.towerTopNodeDof = towerTopNodeDof;
towerModel.towerTopSelectionMatrix = towerTopSelectionMatrix;
towerModel.massMatrixRNA_atTowerTop = massMatrixRNA_atTowerTop;
towerModel.K_global = K_global;
towerModel.M_global = M_global;
towerModel.M_tower_global = M_tower_global;
towerModel.M_RNA_global = M_RNA_global;
end
