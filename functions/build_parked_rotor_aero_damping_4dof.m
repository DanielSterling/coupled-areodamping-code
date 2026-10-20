function rotorOut = build_parked_rotor_aero_damping_4dof( ...
    bladeAeroData, rnaOut, towerTop_xyz, operatingPoint, ...
    aerofoilDatabase, aerofoilSectionSet)
%BUILD_PARKED_ROTOR_AERO_DAMPING_4DOF 停机三叶片塔顶4x4气动阻尼。
%
% 物理对象：柔性塔筒上的刚性RNA；叶片固定方位、固定独立桨距，转速0。
% 全局坐标：X=SS，Y=FA，Z向上。所有力矩均关于塔顶参考点。
% 输入：
%   bladeAeroData：read_blade_aero_geometry的输出，沿用blade_geo.xlsx。
%   rnaOut：build_rna_lumped_mass_at_tower_top的输出，复用叶根位置和旋转。
%   towerTop_xyz：塔顶全局坐标，3x1或1x3，m。
%   operatingPoint：windVelocity_global(3x1, m/s)、airDensity、numGaussPoint。
%   aerofoilDatabase、aerofoilSectionSet：沿用现有气动数据库。
% 输出：
%   rotorOut.dampingMatrix4：三叶片总矩阵，4x4，允许非对称。
%   rotorOut.dampingMatrix4_byBlade：逐叶片矩阵，4x4x3。
%   rotorOut.bladeAeroOut：攻角、系数、截面C2、积分点B4和力臂。
%   rotorOut.bladeSectionTable：逐叶片截面表。
%
% 自由度和广义力严格对应：
%   q4=[uX;uY;thetaX;thetaY]，Q4=[FX;FY;MX;MY]。
%   delta_Q4=-C4*q4dot；装入方程左侧使用+C4*q4dot。
%   thetaX对应FA弯曲转角，thetaY对应SS弯曲转角。
%
% 核心物理过程：
%   v_global=G4*q4dot，G4=[1 0 0 rz;0 1 -rz 0;0 0 ry -rx]。
%   B4=R_reference(:,1:2).'*G4；截面速度vXY=B4*q4dot。
%   Q4=B4.'*fXY；C4=sum_blades integral(B4.'*Ca2*B4 ds)。
%   直接构造B4和C4，不先生成6x6矩阵再截取。
%
% 范围：仅保留四个塔顶气动通道；忽略uZ、thetaZ速度反馈及FZ、MZ载荷。
% 不约束结构的uZ/thetaZ；结构是否保留这些自由度由结构模型另行决定。
% 包括倾角、锥角、方位角、叶根偏置；桨距/扭角只进入攻角。
% 直线叶片，气动力作用于中心线；不含Cm、叶片柔性、平均载荷、
% 气动刚度、非定常气动及塔筒自身涡激载荷。

%% 1. 统一物理量与积分规则
towerTop_xyz = towerTop_xyz(:);
windVelocity_global = operatingPoint.windVelocity_global(:);
numBlade = numel(rnaOut.threeBlades.blade);
numSection = numel(bladeAeroData.chord);
distanceAlongBlade = bladeAeroData.distanceAlongBlade(:);
[gaussCoordinate, gaussWeight] = gauss_rule(operatingPoint.numGaussPoint);
dampingMatrix4_byBlade = zeros(4, 4, numBlade);
bladeSectionTable = cell(numBlade, 1);

for bladeIdx = 1:numBlade
    %% 2. 移除RNA桨距旋转，得到截面公式使用的参考坐标
    geometry = rnaOut.threeBlades.blade(bladeIdx);
    pitch = geometry.pitch_rad;
    rotationPitch = [cos(pitch), sin(pitch), 0; ...
                    -sin(pitch), cos(pitch), 0; 0, 0, 1];
    rotationReferenceToGlobal = geometry.rotRootToGlobal * rotationPitch.';
    referenceAxesXY_global = rotationReferenceToGlobal(:, 1:2);
    windVelocity_reference = rotationReferenceToGlobal.' * windVelocity_global;
    Vx = windVelocity_reference(1);
    Vy = windVelocity_reference(2);

    %% 3. 每支叶片分别查询当前桨距下的极曲线及对弧度导数
    sectionAoA_rad = atan2(Vx, Vy) - pitch - bladeAeroData.aeroTwist_rad(:);
    sectionAoA_rad = mod(sectionAoA_rad + pi, 2*pi) - pi;
    [Cl, Cd, ~, dCl_dalpha, dCd_dalpha, ~] = query_aerofoil_section( ...
        aerofoilDatabase, aerofoilSectionSet, bladeAeroData.foilSection, ...
        bladeAeroData.thickness_pct, sectionAoA_rad);
    sectionDampingMatrix2 = zeros(2, 2, numSection);
    for sectionIdx = 1:numSection
        sectionOut = aero_damping_section_scalar(Vx, Vy, ...
            bladeAeroData.chord(sectionIdx), operatingPoint.airDensity, ...
            Cl(sectionIdx), Cd(sectionIdx), ...
            dCl_dalpha(sectionIdx), dCd_dalpha(sectionIdx));
        sectionDampingMatrix2(:, :, sectionIdx) = sectionOut.Ca2;
    end

    %% 4. 截面中心线位置；RNA叶根位置已经包含轮毂到叶根的偏置
    sectionPosition_global = geometry.rootReference_xyz(:) ...
        + rotationReferenceToGlobal * [zeros(2, numSection); ...
        bladeAeroData.spanPositionFromRoot(:).'];

    %% 5. 直接以2x4速度映射积分单叶片塔顶4x4矩阵
    numPoint = (numSection - 1) * numel(gaussCoordinate);
    gaussPointInfo.segmentIndex = zeros(numPoint, 1);
    gaussPointInfo.distanceAlongBlade = zeros(numPoint, 1);
    gaussPointInfo.integrationWeight_m = zeros(numPoint, 1);
    gaussPointInfo.rTowerTopToSection_global = zeros(3, numPoint);
    gaussPointInfo.sectionVelocityMap4 = zeros(2, 4, numPoint);
    gaussPointInfo.sectionDampingMatrix2 = zeros(2, 2, numPoint);
    gaussPointInfo.dampingContribution4 = zeros(4, 4, numPoint);
    dampingMatrix4 = zeros(4);
    pointIdx = 0;
    for segmentIdx = 1:numSection-1
        segmentLength = distanceAlongBlade(segmentIdx+1) ...
            - distanceAlongBlade(segmentIdx);
        for gaussIdx = 1:numel(gaussCoordinate)
            pointIdx = pointIdx + 1;
            shape1 = (1 - gaussCoordinate(gaussIdx))/2;
            shape2 = (1 + gaussCoordinate(gaussIdx))/2;
            sectionPosition = shape1*sectionPosition_global(:, segmentIdx) ...
                + shape2*sectionPosition_global(:, segmentIdx+1);
            leverArm = sectionPosition - towerTop_xyz;
            sectionMatrix = shape1*sectionDampingMatrix2(:, :, segmentIdx) ...
                + shape2*sectionDampingMatrix2(:, :, segmentIdx+1);

            % 即使不保留塔顶uZ，thetaX/thetaY仍会产生截面竖向速度。
            % 第三行不能删除；倾角/锥角可能将其投影进截面XY平面。
            rx = leverArm(1); ry = leverArm(2); rz = leverArm(3);
            rigidVelocityMap4 = [1, 0, 0, rz; ...
                                 0, 1, -rz, 0; ...
                                 0, 0, ry, -rx];
            sectionVelocityMap4 = referenceAxesXY_global.' * rigidVelocityMap4;
            integrationWeight = segmentLength/2 * gaussWeight(gaussIdx);
            contribution4 = sectionVelocityMap4.' * sectionMatrix ...
                * sectionVelocityMap4 * integrationWeight;
            dampingMatrix4 = dampingMatrix4 + contribution4;

            gaussPointInfo.segmentIndex(pointIdx) = segmentIdx;
            gaussPointInfo.distanceAlongBlade(pointIdx) = ...
                shape1*distanceAlongBlade(segmentIdx) ...
                + shape2*distanceAlongBlade(segmentIdx+1);
            gaussPointInfo.integrationWeight_m(pointIdx) = integrationWeight;
            gaussPointInfo.rTowerTopToSection_global(:, pointIdx) = leverArm;
            gaussPointInfo.sectionVelocityMap4(:, :, pointIdx) = sectionVelocityMap4;
            gaussPointInfo.sectionDampingMatrix2(:, :, pointIdx) = sectionMatrix;
            gaussPointInfo.dampingContribution4(:, :, pointIdx) = contribution4;
        end
    end

    %% 6. 保存单叶片贡献，所有叶片采用相同塔顶参考点与全局轴
    currentBlade = struct();
    currentBlade.bladeIdx = bladeIdx;
    currentBlade.bladeAzimuth_rad = geometry.azimuth_rad;
    currentBlade.bladePitch_rad = pitch;
    currentBlade.dampingMatrix4 = dampingMatrix4;
    currentBlade.rotationReferenceToGlobal = rotationReferenceToGlobal;
    currentBlade.rootReference_xyz = geometry.rootReference_xyz(:);
    currentBlade.sectionPosition_global = sectionPosition_global;
    currentBlade.windVelocity_global = windVelocity_global;
    currentBlade.windVelocity_reference = windVelocity_reference;
    currentBlade.sectionAoA_rad = sectionAoA_rad;
    currentBlade.sectionLiftCoeff = Cl;
    currentBlade.sectionDragCoeff = Cd;
    currentBlade.sectionLiftCoeff_dAlpha = dCl_dalpha;
    currentBlade.sectionDragCoeff_dAlpha = dCd_dalpha;
    currentBlade.sectionDampingMatrix2 = sectionDampingMatrix2;
    currentBlade.gaussPointInfo = gaussPointInfo;
    if bladeIdx == 1
        bladeAeroOut = repmat(currentBlade, numBlade, 1);
    else
        bladeAeroOut(bladeIdx) = currentBlade;
    end
    dampingMatrix4_byBlade(:, :, bladeIdx) = dampingMatrix4;
    bladeSectionTable{bladeIdx} = table(distanceAlongBlade, ...
        bladeAeroData.spanPositionFromRoot(:), bladeAeroData.chord(:), ...
        rad2deg(sectionAoA_rad), Cl, Cd, dCl_dalpha, dCd_dalpha, ...
        'VariableNames', {'DistanceAlongBlade_m', 'RootZ_m', 'Chord_m', ...
        'AoA_deg', 'Cl', 'Cd', 'dCl_dAlpha_perRad', 'dCd_dAlpha_perRad'});
end

%% 7. 三叶片求和；不施加运行风轮的对称抵消条件
rotorOut.dampingMatrix4 = sum(dampingMatrix4_byBlade, 3);
rotorOut.dampingMatrix4_byBlade = dampingMatrix4_byBlade;
rotorOut.bladeAeroOut = bladeAeroOut;
rotorOut.bladeSectionTable = bladeSectionTable;
rotorOut.topDofIdIn6 = [1; 2; 4; 5];
rotorOut.displacementOrder = {'uX_SS', 'uY_FA', 'thetaX', 'thetaY'};
rotorOut.loadOrder = {'FX_SS', 'FY_FA', 'MX', 'MY'};
rotorOut.rotorSpeed_rpm = 0;
rotorOut.operatingPoint = operatingPoint;
rotorOut.towerTop_xyz = towerTop_xyz;
rotorOut.modelDescription = ...
    '停机刚性RNA；仅四个塔顶气动通道，省略FZ/MZ和uZ/thetaZ气动反馈';
end

function [coordinate, weight] = gauss_rule(numPoint)
switch numPoint
    case 2
        coordinate = [-1; 1]/sqrt(3);
        weight = [1; 1];
    case 3
        coordinate = [-sqrt(3/5); 0; sqrt(3/5)];
        weight = [5/9; 8/9; 5/9];
    otherwise
        error('numGaussPoint支持2或3。');
end
end
