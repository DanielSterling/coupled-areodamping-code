%% 停机塔顶四自由度气动阻尼：直接构造4x4，可接续现有结构时程
% 物理对象：柔性塔筒 + 随塔顶刚体运动的机舱、轮毂和三支叶片。
% 全局坐标：X = SS 侧向；Y = FA 前后；Z = 竖直向上。
% 每个塔筒节点自由度：[uX uY uZ thetaX thetaY thetaZ]。
% 功能：结构建模、三叶片塔顶气动阻尼、模态截断、短时扰动力激励及撤载后的时程求解。
% 每支叶片Ca2 -> 2x4速度映射B4 -> 塔顶Ca4 -> 三叶片求和 -> P4装配。
% 四个气动自由度[uX,uY,thetaX,thetaY]；结构仍为每节点6自由度。
% 省略FZ/MZ及uZ/thetaZ气动反馈，不会把结构uZ/thetaZ约束为零。
% RNA 内部不设置柔性自由度；气动阻尼在当前固定工况下线性化。
% 当前仅包含叶片气动阻尼反馈，不计算平均气动力、静平衡或气动刚度。
% 时程输出为所采用参考状态附近的扰动响应；未按静变形更新气动工况。
% 初始位移和速度为0；塔顶施加有限周期的平滑谐波扰动力，随后撤载。
% 从机型数据开始构建模型；扰动力为结构外载，固定均匀风用于计算气动阻尼。
% 沿用现有 build_rna_lumped_mass_at_tower_top 及其数据格式。
clear
clc

%% 1. 项目路径与计算设置
projectFolder = fileparts(mfilename('fullpath'));
addpath(fullfile(projectFolder, 'functions'));
turbineModelName = "my260";
turbineDataFolder = fullfile(projectFolder, turbineModelName);
towerFile = fullfile(turbineDataFolder, "tatong.xlsx");
nacelleFile = fullfile(turbineDataFolder, "nacelle.xlsx");
hubFile = fullfile(turbineDataFolder, "hub.xlsx");
bladeRNAFile = fullfile(turbineDataFolder, "blade_rna.xlsx");
bladeGeoFile = fullfile(turbineDataFolder, "blade_geo.xlsx");
aerofoilDatabaseFile = fullfile(turbineDataFolder, "AerofoilDatabase.xlsx");
aerofoilSectionFile = fullfile(turbineDataFolder, "AerofoilSectionSet.xlsx");
% 停机风轮几何设置；锥角、倾角沿用 RNA 函数从文件读取的值。
bladeAzimuthOffset_rad = 0;
bladeAzimuthSet_rad = bladeAzimuthOffset_rad + deg2rad([0; 120; 240]);
bladePitchSet_rad = deg2rad([90; 90; 90]);
numModesRequested = 16;                   % 求解的无阻尼模态数
numModesRetained = 6;                     % 参与耦合计算的模态数
enableTimeResponse = true;              % false：计算并保存C4后结束；true：继续耦合时程
structuralDampingRatio = 0;               % 标量，或保留模态对应的列向量
% 塔顶扰动力设置；六分量顺序[FX;FY;FZ;MX;MY;MZ]，单位[N;N;N;N*m;N*m;N*m]。
excitationAmplitude6 = [0; 2.0e4; 0; 0; 0; 0]; % 默认仅施加FA方向扰动力
excitationStartTime_s = 0;              % 开始施力的时刻，s；对齐输出网格
excitationNumCycles = 4;                 % 施力周期数
excitationFrequency_Hz = [];             % []：自动选主要FA模态频率；也可填写0.25
% 载荷形式：amplitude6 * sin(2*pi*f*tau) * sin(pi*tau/T)^2。
% amplitude6为正弦载荷的幅值参数；乘窗后采样峰值通常略小于该值。
timeStep_s = 0.02;                        % 时程输出间隔，s
totalTime_s = 400;                        % 时程总时长，s
enableModePlot = false;
modeIndexToPlot = [1; 2];
% 本阶段只计算一个确定工况，三支叶片分别使用各自的方位角和桨距。
inflowSpeed = 6;                          % m/s
inflowAngle_rad = deg2rad(240);           % 空气速度方向，沿用原工况
% inflowAngle_rad = deg2rad(90);  
airDensity = 1.29;                       % kg/m^3
numGaussPointBlade = 3;
% 桨距直接使用上面的bladePitchSet_rad，不设置第二套桨距参数。
% 例如独立变桨：bladePitchSet_rad = deg2rad([78; 90; 85]);
resultFolder = fullfile(projectFolder, "results", "parked_aero_4dof");
resultFile = fullfile(resultFolder, "three_blade_aero_damping_4dof.mat");
displacementFile = fullfile(resultFolder, "tower_top_force_displacement.txt");
responseFile = fullfile(resultFolder, "tower_rna_force_response.mat");
forceFile = fullfile(resultFolder, "tower_top_disturbance_force.txt");

%% 2. 读取塔筒数据并建立梁单元
[towerRawData, ~] = read_xlsx_numeric(towerFile, 1);
[towerNodeData, towerElemData, towerConn] = ...
    tower_elem_props_from_Tatong_gauss_FE(towerRawData);
numTowerElem = numel(towerElemData);
numTowerNode = size(towerNodeData.xyz, 1);
% 沿用塔筒前处理的节点顺序：底部在前，塔顶为最后一个节点。
towerTopNodeId = numTowerNode;
towerTop_xyz = towerNodeData.xyz(towerTopNodeId, :);

%% 3. 构建 RNA 在塔顶、全局坐标系下的 6 x 6 集中质量矩阵
% RNA 函数内部已包含三个部件的质量、转动惯量及偏置耦合。
% 此处只调用一次，不再单独叠加机舱、轮毂或叶片质量。
rnaOut = build_rna_lumped_mass_at_tower_top( ...
    nacelleFile, hubFile, bladeRNAFile, ...
    towerTop_xyz, bladeAzimuthSet_rad, bladePitchSet_rad);
massMatrixRNA_atTowerTop = rnaOut.massMatrix6;

%% 4. 装配塔筒 + RNA 总体质量矩阵和刚度矩阵
% M_global = M_tower_global + P' * massMatrixRNA_atTowerTop * P
% K_global = K_tower_global
% P 从全部塔筒自由度中提取塔顶的六个自由度。
towerModel = assemble_tower_rna_structure( ...
    towerNodeData, towerElemData, towerConn, massMatrixRNA_atTowerTop);
K_global = towerModel.K_global;
M_global = towerModel.M_global;
M_tower_global = towerModel.M_tower_global;
M_RNA_global = towerModel.M_RNA_global;
dofPerNode = towerModel.dofPerNode;
numDofTotal = towerModel.numDofTotal;
elemDofMap = towerModel.elemDofMap;
towerTopNodeDof = towerModel.towerTopNodeDof;

%% 5. 塔底六自由度固定，提取自由自由度子矩阵
fixedDofId = (1:dofPerNode).';
freeDofId = setdiff((1:numDofTotal).', fixedDofId);
K_free = K_global(freeDofId, freeDofId);
M_free = M_global(freeDofId, freeDofId);
% 后续的塔顶气动力、气动阻尼和响应提取可直接使用这个映射：
%   q_top = towerTopSelectionMatrix_free * q_free
%   F_free = towerTopSelectionMatrix_free' * F_top
%   C_free = towerTopSelectionMatrix_free' * C_top * ...
%            towerTopSelectionMatrix_free
towerTopSelectionMatrix_free = ...
    towerModel.towerTopSelectionMatrix(:, freeDofId);
[~, towerTopFreeDofId] = ismember(towerTopNodeDof, freeDofId);

%% 6. 求解无阻尼固有频率和质量归一化振型
% K_free * Phi = M_free * Phi * diag(omega.^2)
% Phi' * M_free * Phi = I
modalOut = solve_tower_rna_undamped_modes( ...
    K_free, M_free, numModesRequested);
modeShapeUndamped = modalOut.modeShapeUndamped;
naturalOmegaUndamped_rad_s = modalOut.naturalOmegaUndamped_rad_s;
naturalFreqUndamped_Hz = modalOut.naturalFreqUndamped_Hz;
% 恢复塔底零位移，得到完整塔筒节点振型。
modeShapeUndamped_global = zeros(numDofTotal, modalOut.numModesSolved);
modeShapeUndamped_global(freeDofId, :) = modeShapeUndamped;
towerTopModeShape = towerTopSelectionMatrix_free * modeShapeUndamped;

%% 7. 整理结构模型输出并显示频率
% dynamicOut 是本阶段结果容器，主程序的各个物理步骤仍直接可见。
dynamicOut = towerModel;
dynamicOut.towerNodeData = towerNodeData;
dynamicOut.towerElemData = towerElemData;
dynamicOut.towerTop_xyz = towerTop_xyz;
dynamicOut.rnaOut = rnaOut;
dynamicOut.bladeAzimuthSet_rad = bladeAzimuthSet_rad;
dynamicOut.bladePitchSet_rad = bladePitchSet_rad;
dynamicOut.fixedDofId = fixedDofId;
dynamicOut.freeDofId = freeDofId;
dynamicOut.towerTopFreeDofId = towerTopFreeDofId;
dynamicOut.towerTopSelectionMatrix_free = towerTopSelectionMatrix_free;
dynamicOut.K_free = K_free;
dynamicOut.M_free = M_free;
dynamicOut.modeShapeUndamped = modeShapeUndamped;
dynamicOut.modeShapeUndamped_global = modeShapeUndamped_global;
dynamicOut.towerTopModeShape = towerTopModeShape;
dynamicOut.naturalOmegaUndamped_rad_s = naturalOmegaUndamped_rad_s;
dynamicOut.naturalFreqUndamped_Hz = naturalFreqUndamped_Hz;
dynamicOut.numModesRequested = numModesRequested;
dynamicOut.numModesSolved = modalOut.numModesSolved;
dynamicOut.eigSolveFlag = modalOut.eigSolveFlag;
modeFrequencyTable = table( ...
    (1:modalOut.numModesSolved).', naturalFreqUndamped_Hz, ...
    'VariableNames', {'ModeIndex', 'Frequency_Hz'});
fprintf('塔筒单元数：%d；塔筒节点数：%d\n', numTowerElem, numTowerNode);
fprintf('总体自由度：%d；约束后自由度：%d\n', ...
    numDofTotal, numel(freeDofId));
disp('RNA 在塔顶的 6 x 6 质量矩阵：');
disp(massMatrixRNA_atTowerTop);
disp('塔筒 + RNA 无阻尼固有频率：');
disp(modeFrequencyTable);

%% 8. 可选：查看塔筒三个平动方向的振型分量
% 振型是相对运动形状，不是真实响应幅值。
if enableModePlot
    plot_tower_rna_modes(towerNodeData, ...
        modeShapeUndamped_global, naturalFreqUndamped_Hz, modeIndexToPlot);
end

%% 9. 读取翼型数据库和叶片气动截面数据
% 气动截面读取独立于叶片柔性建模，不需要blade_mass_stiff.xlsx。
aerofoilDatabase = build_aerofoil_database(aerofoilDatabaseFile);
aerofoilSectionSet = build_aerofoil_section_set( ...
    aerofoilDatabase, aerofoilSectionFile);
bladeAeroData = read_blade_aero_geometry(bladeGeoFile);

%% 10. 构建三支叶片共用的停机空气速度工况
% 均匀空气速度相同；各叶片在自己的参考坐标系内重新投影。
windVelocity_global = inflowSpeed ...
    * [cos(inflowAngle_rad); sin(inflowAngle_rad); 0];
operatingPoint.windVelocity_global = windVelocity_global;
operatingPoint.inflowSpeed = inflowSpeed;
operatingPoint.inflowAngle_rad = inflowAngle_rad;
operatingPoint.airDensity = airDensity;
operatingPoint.numGaussPoint = numGaussPointBlade;
operatingPoint.bladeAzimuthSet_rad = bladeAzimuthSet_rad;
operatingPoint.bladePitchSet_rad = bladePitchSet_rad;

%% 11. 直接计算停机风轮塔顶四自由度气动阻尼
% q4=[uX;uY;thetaX;thetaY]，Q4=[FX;FY;MX;MY]。
% 不先计算6x6；直接积分C4=sum integral(B4.'*Ca2*B4 ds)。
% 倾角、锥角、方位角、独立桨距与叶根位置沿用RNA和气动数据。
rotorAero4Out = build_parked_rotor_aero_damping_4dof( ...
    bladeAeroData, rnaOut, towerTop_xyz, operatingPoint, ...
    aerofoilDatabase, aerofoilSectionSet);
numBlade = numel(rotorAero4Out.bladeAeroOut);
bladeAeroOut = rotorAero4Out.bladeAeroOut;
bladeSectionTable = rotorAero4Out.bladeSectionTable;
dampingMatrixBlade4_atTowerTop = rotorAero4Out.dampingMatrix4_byBlade;
dampingMatrixRotor4_atTowerTop = rotorAero4Out.dampingMatrix4;
topAeroDofIdIn6 = rotorAero4Out.topDofIdIn6;

%% 12. 四通道矩阵的方向与单位
% 列：uXdot,uYdot,thetaXdot,thetaYdot；行：FX,FY,MX,MY。
% 左上2x2单位N*s/m；交叉块N*s；右下2x2单位N*m*s（rad无量纲）。
% 无旋转风轮的三叶片一般不能使用运行工况的零元素分布；保留全部16项。

%% 13. 将4x4矩阵直接装配到现有六自由度梁结构
% P4取塔顶全局[uX,uY,thetaX,thetaY]；C_global=P4.'*C4*P4。
% 这是四通道气动近似；结构仍允许uZ/thetaZ运动，但不计算其气动反馈。
% 不进行静力凝聚，不用Schur补，也不把结构自由度数改成4。
towerTopAeroSelectionMatrix = ...
    towerModel.towerTopSelectionMatrix(topAeroDofIdIn6, :);
towerTopAeroSelectionMatrix_free = towerTopAeroSelectionMatrix(:, freeDofId);
C_aero_tower_global = towerTopAeroSelectionMatrix.' ...
    * sparse(dampingMatrixRotor4_atTowerTop) * towerTopAeroSelectionMatrix;
C_aero_tower_free = C_aero_tower_global(freeDofId, freeDofId);

%% 14. 显示逐叶片与总4x4矩阵
% dampingMatrixRotor4_atTowerTop是本次需要的物理空间塔顶4x4矩阵。
dynamicOut.operatingPoint = operatingPoint;
dynamicOut.rotorAero4Out = rotorAero4Out;
dynamicOut.bladeAeroOut = bladeAeroOut;
dynamicOut.dampingMatrixBlade4_atTowerTop = dampingMatrixBlade4_atTowerTop;
dynamicOut.dampingMatrixRotor4_atTowerTop = dampingMatrixRotor4_atTowerTop;
dynamicOut.towerTopAeroSelectionMatrix_free = towerTopAeroSelectionMatrix_free;
dynamicOut.C_aero_tower_global = C_aero_tower_global;
dynamicOut.C_aero_tower_free = C_aero_tower_free;
fprintf('停机0 rpm；空气速度%.2f m/s；方向角%.2f deg\n', ...
    inflowSpeed, rad2deg(inflowAngle_rad));
for bladeIdx = 1:numBlade
    fprintf('叶片%d：方位%.2f deg；桨距%.2f deg\n', bladeIdx, ...
        rad2deg(bladeAeroOut(bladeIdx).bladeAzimuth_rad), ...
        rad2deg(bladeAeroOut(bladeIdx).bladePitch_rad));
    disp('单叶片C4：行[FX,FY,MX,MY]，列[uXdot,uYdot,thetaXdot,thetaYdot]');
    disp(dampingMatrixBlade4_atTowerTop(:, :, bladeIdx));
end
topAeroDofLabels = {'uX_SS', 'uY_FA', 'thetaX', 'thetaY'};
topAeroLoadLabels = {'FX_SS', 'FY_FA', 'MX', 'MY'};
dampingMatrix4Table = array2table(dampingMatrixRotor4_atTowerTop, ...
    'VariableNames', topAeroDofLabels, 'RowNames', topAeroLoadLabels);
disp('三叶片风轮总塔顶4x4气动阻尼矩阵：');
disp(dampingMatrix4Table);

%% 15. 保存三支叶片矩阵、风轮总矩阵及塔筒气动阻尼矩阵
% MAT文件保存稀疏矩阵及完整工况信息；同名文件在下次运行时更新。
% 同时保存结构模型和模态；本主程序随后继续完成模态截断与扰动力时程。
if ~isfolder(resultFolder)
    mkdir(resultFolder);
end
save(resultFile, ...
    'rotorAero4Out', 'topAeroDofIdIn6', 'dampingMatrix4Table', ...
    'towerTopAeroSelectionMatrix_free', ...
    'dampingMatrixBlade4_atTowerTop', 'dampingMatrixRotor4_atTowerTop', ...
    'C_aero_tower_global', 'C_aero_tower_free', ...
    'bladeAeroOut', 'bladeSectionTable', 'bladeAeroData', 'operatingPoint', ...
    'bladeAzimuthSet_rad', 'bladePitchSet_rad', ...
    'towerTop_xyz', 'towerTopNodeDof', 'towerTopFreeDofId', ...
    'freeDofId', 'fixedDofId', 'towerTopSelectionMatrix_free', ...
    'dynamicOut', 'turbineModelName', '-v7.3');
fprintf('三叶片气动阻尼结果已保存：%s\n', char(resultFile));

% CSV含行列物理量名称，便于检查与对比。
writetable(dampingMatrix4Table, ...
    fullfile(resultFolder, 'rotor_tower_top_aero_damping_4dof.csv'), ...
    'WriteRowNames', true);
if ~enableTimeResponse
    fprintf('4x4气动阻尼计算完成；enableTimeResponse=true可继续时程。\n');
    return
end

%% 步骤5.1 选择保留的无阻尼模态
% 参数统一在程序开头设置；阻尼投影保留非对角项及不对称性。
retainedModeId = (1:numModesRetained).';

%% 步骤5.2 投影质量、刚度和完整气动阻尼矩阵
modeShapeRetained = modeShapeUndamped(:, retainedModeId);
M_modal = modeShapeRetained.' * M_free * modeShapeRetained;
K_modal = modeShapeRetained.' * K_free * modeShapeRetained;
C_aero_modal = modeShapeRetained.' * C_aero_tower_free * modeShapeRetained;

% 无阻尼模态已经质量归一化；结构阻尼按保留模态分别指定。
naturalOmegaRetained_rad_s = naturalOmegaUndamped_rad_s(retainedModeId);
C_struct_modal = diag(2 * structuralDampingRatio(:) ...
    .* naturalOmegaRetained_rad_s .* diag(M_modal));
C_total_modal = C_struct_modal + C_aero_modal;

%% 步骤5.3 构建载荷映射、状态方程和响应恢复映射
% eta为模态位移；q_top=H*eta；F_modal=H.'*F_top。
towerTopModeShapeRetained = towerTopSelectionMatrix_free * modeShapeRetained;
towerTopAeroModeShapeRetained = towerTopAeroSelectionMatrix_free * modeShapeRetained;
forceMatrixTop_modal = towerTopModeShapeRetained.';
A_modal = [zeros(numModesRetained), eye(numModesRetained); ...
    -(M_modal \ K_modal), -(M_modal \ C_total_modal)];
B_top_modal = [zeros(numModesRetained, 6); M_modal \ forceMatrixTop_modal];
responseMatrixTop_modal = [towerTopModeShapeRetained, zeros(6, numModesRetained)];

modalModelOut = struct();
modalModelOut.numModesRetained = numModesRetained;
modalModelOut.retainedModeId = retainedModeId;
modalModelOut.modeShapeRetained = modeShapeRetained;
modalModelOut.massMatrixModal = M_modal;
modalModelOut.stiffnessMatrixModal = K_modal;
modalModelOut.aeroDampingMatrixModal = C_aero_modal;
modalModelOut.structuralDampingMatrixModal = C_struct_modal;
modalModelOut.structuralDampingRatio = structuralDampingRatio;
modalModelOut.dampingMatrixModal = C_total_modal;
modalModelOut.naturalOmegaRetained_rad_s = naturalOmegaRetained_rad_s;
modalModelOut.naturalFrequencyRetained_Hz = naturalOmegaRetained_rad_s / (2*pi);
modalModelOut.towerTopModeShape = towerTopModeShapeRetained;
modalModelOut.forceMatrixTop_modal = forceMatrixTop_modal;
modalModelOut.stateMatrix = A_modal;
modalModelOut.inputMatrixTop = B_top_modal;
modalModelOut.outputMatrixTop = responseMatrixTop_modal;

%% 步骤5.4 保存逐支叶片的完整模态阻尼贡献
% 通过同一套塔顶振型投影，保留各叶片对所有模态耦合项的贡献。
numBlade = size(dampingMatrixBlade4_atTowerTop, 3);
C_aero_modal_byBlade = zeros(numModesRetained, numModesRetained, numBlade);
for bladeIdx = 1:numBlade
    C_aero_modal_byBlade(:, :, bladeIdx) = ...
        towerTopAeroModeShapeRetained.' ...
        * dampingMatrixBlade4_atTowerTop(:, :, bladeIdx) ...
        * towerTopAeroModeShapeRetained;
end
modalModelOut.aeroDampingMatrixModal_byBlade = C_aero_modal_byBlade;
modalModelOut.topAeroDofIdIn6 = topAeroDofIdIn6;
modalModelOut.towerTopAeroModeShape = towerTopAeroModeShapeRetained;
modalModelOut.rotorAero4Out = rotorAero4Out;

%% 步骤5.5 补齐工况信息并保存耦合模型
% 物理空间塔顶矩阵与模态空间矩阵分别保存，便于后续读取复用。
modalModelOut.operatingPoint = operatingPoint;
modalModelOut.turbineModelName = turbineModelName;
modalModelOut.rnaOut = rnaOut;
modalModelOut.bladeAzimuthSet_rad = bladeAzimuthSet_rad;
modalModelOut.bladePitchSet_rad = bladePitchSet_rad;
modalModelOut.towerTop_xyz = towerTop_xyz;
modalModelOut.freeDofId = freeDofId;
modalModelOut.fixedDofId = fixedDofId;
modalModelOut.towerTopNodeDof = towerTopNodeDof;
modalModelOut.towerTopFreeDofId = towerTopFreeDofId;
modalModelOut.towerTopSelectionMatrix_free = towerTopSelectionMatrix_free;
modalModelOut.dampingMatrixBlade4_atTowerTop = dampingMatrixBlade4_atTowerTop;
modalModelOut.dampingMatrixRotor4_atTowerTop = dampingMatrixRotor4_atTowerTop;
modalModelOut.coordinateDescription = 'X=SS; Y=FA; Z=up';
modalModelOut.referenceStateDescription = ...
    '固定工况下的线性扰动模型；未计算平均气动力、静平衡及气动刚度';
% 模型状态[eta;eta_dot]；塔顶输出顺序[uX,uY,uZ,thetaX,thetaY,thetaZ]。
% 塔顶输入顺序[FX,FY,FZ,MX,MY,MZ]。
% 本段保存模型；时程在步骤6求解。当前尚未加入塔筒自身气动力或气动刚度。
step5ResultFolder = resultFolder;
if ~isfolder(step5ResultFolder)
    mkdir(step5ResultFolder);
end
save(fullfile(step5ResultFolder, 'tower_rna_modal_model.mat'), ...
    'modalModelOut', 'retainedModeId', 'structuralDampingRatio', ...
    'modeShapeRetained', 'M_modal', 'K_modal', ...
    'C_aero_modal', 'C_struct_modal', 'C_total_modal', ...
    'C_aero_modal_byBlade', 'towerTopModeShapeRetained', ...
    'towerTopAeroModeShapeRetained', 'topAeroDofIdIn6', ...
    'forceMatrixTop_modal', 'A_modal', 'B_top_modal', ...
    'responseMatrixTop_modal');
fprintf('步骤5耦合模型已保存：%s\n', ...
    fullfile(step5ResultFolder, 'tower_rna_modal_model.mat'));

%% 步骤6.1 选择扰动力参考频率并构造塔顶载荷
% eta为模态位移，塔顶载荷F_top映射为F_modal=H.'*F_top。
% 自动选频使用各保留模态的FA静柔度贡献H(2,j)^2/K_modal(j,j)。
% 这里只选择激励频率；求解中仍保留全部模态和完整阻尼耦合项。
time_s = (0:timeStep_s:totalTime_s).';
if isempty(excitationFrequency_Hz)
    faFlexibilityContribution = towerTopModeShapeRetained(2, :).'.^2 ...
        ./ diag(K_modal);
    [~, faLocalModeIdx] = max(faFlexibilityContribution);
    targetFrequency_Hz = naturalOmegaRetained_rad_s(faLocalModeIdx) / (2*pi);
else
    targetFrequency_Hz = excitationFrequency_Hz;
end

% 将施力起点和持续时间对齐时间网格，确保完整施加指定周期数。
% 实际频率会相对目标值有轻微调整，二者均保存并打印。
excitation.amplitude6 = excitationAmplitude6;
excitation.startTime_s = round(excitationStartTime_s/timeStep_s)*timeStep_s;
excitation.numCycles = excitationNumCycles;
excitation.targetFrequency_Hz = targetFrequency_Hz;
excitation.duration_s = round(excitationNumCycles ...
    / targetFrequency_Hz / timeStep_s)*timeStep_s;
excitation.frequency_Hz = excitationNumCycles / excitation.duration_s;
excitation.endTime_s = excitation.startTime_s + excitation.duration_s;
excitation.windowDescription = 'sin(pi*tau/T)^2';
if excitation.startTime_s < 0 || excitation.endTime_s > time_s(end)
    error('扰动力时段必须位于仿真时间内，请修改施力起点或总时长。');
end

forceTopHistory = build_tower_top_disturbance(time_s, excitation);
modalForceHistory = forceTopHistory * forceMatrixTop_modal.';

%% 步骤6.2 设置零初始状态，求解完整耦合受迫响应
% M_modal*eta_ddot+C_total_modal*eta_dot+K_modal*eta=H.'*F_top(t)。
% 状态x=[eta;eta_dot]，xdot=A_modal*x+B_top_modal*F_top(t)。
% 气动反馈已包含在C_aero_modal中，整个时程不重新计算气动力。
initialModalDisplacement = zeros(numModesRetained, 1);
initialModalVelocity = zeros(numModesRetained, 1);
initialDisplacementFree = modeShapeRetained * initialModalDisplacement;

% 每个时间步内线性插值外载荷。扩展状态[x;F_top;F_top_dot]后，
% 用矩阵指数精确推进该步内的线性插值载荷；无需求逆A_modal。
numState = 2*numModesRetained;
augmentedMatrix = [A_modal, B_top_modal, zeros(numState, 6); ...
    zeros(6, numState), zeros(6), eye(6); ...
    zeros(6, numState+12)];
transitionMatrix = expm(augmentedMatrix*timeStep_s);
stateTransition = transitionMatrix(1:numState, 1:numState);
forceTransition = transitionMatrix(1:numState, numState+(1:6));
forceSlopeTransition = transitionMatrix(1:numState, numState+6+(1:6));

numTime = numel(time_s);
stateHistory = zeros(numTime, numState);
stateHistory(1, :) = [initialModalDisplacement; initialModalVelocity].';
for timeIdx = 1:numTime-1
    currentForce = forceTopHistory(timeIdx, :).';
    forceSlope = (forceTopHistory(timeIdx+1, :) ...
        - forceTopHistory(timeIdx, :)).'/timeStep_s;
    stateHistory(timeIdx+1, :) = (stateTransition ...
        * stateHistory(timeIdx, :).' + forceTransition*currentForce ...
        + forceSlopeTransition*forceSlope).';
end
modalDisplacement = stateHistory(:, 1:numModesRetained);
modalVelocity = stateHistory(:, numModesRetained+1:end);

%% 步骤6.3 恢复塔顶与塔筒自由自由度响应
% 每行对应一个时刻；塔顶列顺序[uX,uY,uZ,thetaX,thetaY,thetaZ]。
% 计算加速度时包括外部扰动力，撤载后自动退化为自由振动方程。
towerTopDisplacement = modalDisplacement * towerTopModeShapeRetained.';
displacementFA_m = towerTopDisplacement(:, 2);
displacementSS_m = towerTopDisplacement(:, 1);
displacementFree = modalDisplacement * modeShapeRetained.';
towerTopVelocity = modalVelocity * towerTopModeShapeRetained.';
modalAcceleration = (M_modal \ (modalForceHistory.' ...
    - C_total_modal*modalVelocity.' - K_modal*modalDisplacement.')).';
towerTopAcceleration = modalAcceleration * towerTopModeShapeRetained.';

timeResponseOut.time_s = time_s;
timeResponseOut.excitation = excitation;
timeResponseOut.excitationEndTime_s = excitation.endTime_s;
timeResponseOut.forceTopHistory = forceTopHistory;
timeResponseOut.modalForceHistory = modalForceHistory;
timeResponseOut.initialModalDisplacement = initialModalDisplacement;
timeResponseOut.initialModalVelocity = initialModalVelocity;
timeResponseOut.initialDisplacementFree = initialDisplacementFree;
timeResponseOut.stateHistory = stateHistory;
timeResponseOut.modalDisplacement = modalDisplacement;
timeResponseOut.modalVelocity = modalVelocity;
timeResponseOut.modalAcceleration = modalAcceleration;
timeResponseOut.displacementFree = displacementFree;
timeResponseOut.towerTopDisplacement = towerTopDisplacement;
timeResponseOut.towerTopVelocity = towerTopVelocity;
timeResponseOut.towerTopAcceleration = towerTopAcceleration;
timeResponseOut.displacementFA_m = displacementFA_m;
timeResponseOut.displacementSS_m = displacementSS_m;
timeResponseOut.freeDofId = freeDofId;
timeResponseOut.operatingPoint = operatingPoint;

%% 步骤6.4 保存位移TXT、扰动力TXT和完整响应MAT
% 位移TXT：时间(s)、FA位移(m)、SS位移(m)，与初始位移版保持相同列顺序。
% 扰动力TXT：时间(s)、FX、FY、FZ(N)、MX、MY、MZ(N*m)。
% 本TXT为通用数据表，不是Bladed Impact Data File专用格式。
fileId = fopen(displacementFile, 'w');
fprintf(fileId, 'Time_s\tDisplacementFA_m\tDisplacementSS_m\n');
fprintf(fileId, '%.6f\t%.12e\t%.12e\n', ...
    [time_s, displacementFA_m, displacementSS_m].');
fclose(fileId);
fileId = fopen(forceFile, 'w');
fprintf(fileId, 'Time_s\tFX_N\tFY_N\tFZ_N\tMX_Nm\tMY_Nm\tMZ_Nm\n');
fprintf(fileId, '%.6f\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\n', ...
    [time_s, forceTopHistory].');
fclose(fileId);
save(responseFile, 'timeResponseOut', 'modalModelOut', 'excitation', ...
    'timeStep_s', 'totalTime_s', '-v7.3');

%% 步骤6.5 绘制扰动力和塔顶FA、SS位移
responseFigure = figure('Color', 'w', 'Name', '塔顶扰动力响应');
tiledlayout(3, 1);
nexttile;
plot(time_s, forceTopHistory(:, 1:3)/1e3, 'LineWidth', 1.1);
grid on;
xlabel('时间 / s'); ylabel('扰动力 / kN');
legend('FX (SS)', 'FY (FA)', 'FZ');
title('塔顶外部扰动力');
xline(excitation.endTime_s, '--', '撤载');
nexttile;
plot(time_s, displacementFA_m, 'LineWidth', 1.1);
grid on;
xlabel('时间 / s'); ylabel('FA位移 / m');
title('塔顶FA扰动位移');
xline(excitation.endTime_s, '--', '撤载');
nexttile;
plot(time_s, displacementSS_m, 'LineWidth', 1.1);
grid on;
xlabel('时间 / s'); ylabel('SS位移 / m');
title('塔顶SS扰动位移');
xline(excitation.endTime_s, '--', '撤载');
exportgraphics(responseFigure, fullfile(resultFolder, ...
    'tower_top_force_response.png'), 'Resolution', 180);

fprintf('初始位移、初始速度均为0。\n');
fprintf('扰动力时段：%.4f ~ %.4f s；周期数：%g。\n', ...
    excitation.startTime_s, excitation.endTime_s, excitation.numCycles);
fprintf('目标频率：%.8f Hz；实际施力频率：%.8f Hz。\n', ...
    excitation.targetFrequency_Hz, excitation.frequency_Hz);
fprintf('施力期间FA扰动力采样峰值：%.6f kN。\n', ...
    max(abs(forceTopHistory(:, 2)))/1e3);
fprintf('塔顶位移TXT已保存：%s\n', displacementFile);
fprintf('扰动力TXT已保存：%s\n', forceFile);
fprintf('完整响应MAT已保存：%s\n', responseFile);

%% 局部函数：有限周期、sin^2包络的塔顶扰动力
function forceTopHistory = build_tower_top_disturbance(time_s, excitation)
% 输入：时间列向量、载荷起点/持续时间/频率/六分量幅值参数。
% 输出：numTime x 6的塔顶载荷；作用点为塔顶，坐标为全局坐标。
% tau为从施力开始计时的局部时间；施力前后载荷严格为0。
% 平滑包络避免起止时刻发生载荷突跳。
    forceTopHistory = zeros(numel(time_s), 6);
    activeTimeId = time_s > excitation.startTime_s ...
        & time_s < excitation.endTime_s;
    localTime_s = time_s(activeTimeId) - excitation.startTime_s;
    window = sin(pi*localTime_s/excitation.duration_s).^2;
    forceShape = window .* sin(2*pi*excitation.frequency_Hz*localTime_s);
    forceTopHistory(activeTimeId, :) = forceShape * excitation.amplitude6.';
end
