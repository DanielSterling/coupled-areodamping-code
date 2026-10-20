function rnaOut = build_rna_lumped_mass_at_tower_top( ...
    nacelleFile, hubFile, bladeMassFile, towerTop_xyz, ...
    bladeAzimuthSet_rad, bladePitchSet_rad)
%BUILD_RNA_LUMPED_MASS_AT_TOWER_TOP 独立构造完整RNA的塔顶集中质量模型。
%
% 功能：读取机舱、轮毂和叶片的Excel数据，计算各部件质心与惯量，
% 等效到同一塔顶参考点A，叠加得到完整RNA的6×6质量矩阵。
% 所有计算及辅助函数都包含在本文件中，可独立使用。
%
% 输入：
%   nacelleFile：nacelle.xlsx路径，第一列名称、第二列数值。
%     mass [kg]；side、front、top [m]：质心相对塔顶的X、Y、Z偏置；
%     Yaw [kg*m^2]：绕塔筒偏航轴的惯量；
%     Noding、Rolling [kg*m^2]：质心处绕全局X、Y轴的惯量。
%   hubFile：hub.xlsx路径，第一列名称、第二列数值或文本。
%     mass [kg]；centre_of_mass [m]：从轮毂参考点H沿主轴远离塔筒
%     方向到轮毂质心的有符号距离；inertia_shaft、
%     inertia_perpendicular_shaft [kg*m^2]：轮毂质心处轴向和横向惯量；
%     rotor_position [Upwind/Downwind]；tilt_angle_deg [deg]；
%     overhang [m]：H到塔筒中心线的水平距离；
%     hub_vertical_offset、lateral_offset [m]：H相对塔顶的Z、X偏置；
%     cone_angle [deg]：Rotor锥角，正值向下风偏转；
%     blade_root_length [m]：H到叶根沿变桨轴的偏置L。
%   bladeMassFile：blade_rna.xlsx路径，第一列名称、第二列数值。
%     mass [kg]：单支叶片质量；first_mass_moment [kg*m]：叶根一阶
%     质量矩；second_mass_moment [kg*m^2]：叶根二阶质量矩。
%     一、二阶质量矩采用用户指定的同一套叶根长度标准。
%   towerTop_xyz：塔顶全局坐标，3×1或1×3，m。
%   bladeAzimuthSet_rad：三支叶片方位角，三个元素，rad。
%   bladePitchSet_rad：三支叶片桨距角，三个元素，rad。
%
% 坐标与假设：
%   全局X侧向、Y前后且正向指向上风、Z竖直向上。
%   自由度顺序：[Ux, Uy, Uz, Rx, Ry, Rz]。
%   停机、刚性RNA、固定桨距、零机舱偏航；塔顶与偏航轴承参考点重合。
%   机舱质心惯性主轴沿全局坐标轴；轮毂质心惯量绕主轴轴对称。
%   三支叶片采用相同的直线质量分布，质心在展向轴上，忽略预弯、
%   扫掠及截面尺寸产生的自身惯量；Blade mounting cone angle取0。
%   L对应连接段的质量已归入轮毂，不再增加叶片质量。
%   Upwind、零倾角/锥角/方位角/桨距时：x_root=-Y、y_root=X、z_root=Z。
%   方位角绕未锥变的+x_root，桨距绕锥变后+z_root旋转-pitch。
%   当前叶片质心和惯量绕展向轴对称，因此桨距不改变质量矩阵。
%
% 主要输出：
%   massMatrix6：完整RNA关于塔顶的6×6质量矩阵。
%   mass：RNA总质量，kg；massCenter_xyz：总质心全局坐标，3×1，m。
%   rTowerTopToMassCenter：塔顶到总质心的全局偏置，3×1，m。
%   inertiaAtTowerTop_global、inertiaAtMassCenter_global：RNA关于塔顶、
%     总质心的全局惯量张量，3×3，kg*m^2。
%   nacelle、hub、threeBlades：机舱、轮毂和三叶片总结果。
%   threeBlades.blade(k)：第k支叶片的位置、惯量和塔顶质量矩阵。

%% -------------------- 1. 读取三个部件的数据文件 --------------------
nacelleData = readcell(nacelleFile);
nacelleItemName = lower(strtrim(string(nacelleData(:,1))));
nacelleItemValue = str2double(string(nacelleData(:,2)));
nacelleMass = nacelleItemValue(nacelleItemName == "mass");
sideOffset = nacelleItemValue(nacelleItemName == "side");
frontOffset = nacelleItemValue(nacelleItemName == "front");
heightOffset = nacelleItemValue(nacelleItemName == "top");
yawInertia_aboutTowerAxis = nacelleItemValue(nacelleItemName == "yaw");
nodingInertia_aboutMassCenter = nacelleItemValue(nacelleItemName == "noding");
rollingInertia_aboutMassCenter = nacelleItemValue(nacelleItemName == "rolling");

hubData = readcell(hubFile);
hubItemName = lower(strtrim(string(hubData(:,1))));
hubItemText = strtrim(string(hubData(:,2)));
hubMass = str2double(hubItemText(hubItemName == "mass"));
massCenterOffsetAlongShaft = ...
    str2double(hubItemText(hubItemName == "centre_of_mass"));
inertiaShaft = str2double(hubItemText(hubItemName == "inertia_shaft"));
inertiaPerpendicularShaft = ...
    str2double(hubItemText(hubItemName == "inertia_perpendicular_shaft"));
rotorPosition = lower(hubItemText(hubItemName == "rotor_position"));
tiltAngle_rad = deg2rad(str2double(hubItemText(hubItemName == "tilt_angle_deg")));
coneAngle_rad = deg2rad(str2double(hubItemText(hubItemName == "cone_angle")));
overhang = str2double(hubItemText(hubItemName == "overhang"));
hubVerticalOffset = str2double(hubItemText(hubItemName == "hub_vertical_offset"));
lateralOffset = str2double(hubItemText(hubItemName == "lateral_offset"));
bladeRootLength = str2double(hubItemText(hubItemName == "blade_root_length"));

bladeData = readcell(bladeMassFile);
bladeItemName = lower(strtrim(string(bladeData(:,1))));
bladeItemValue = str2double(string(bladeData(:,2)));
bladeMass = bladeItemValue(bladeItemName == "mass");
bladeFirstMassMoment = bladeItemValue(bladeItemName == "first_mass_moment");
bladeSecondMassMoment = bladeItemValue(bladeItemName == "second_mass_moment");

towerTop_xyz = towerTop_xyz(:);
bladeAzimuthSet_rad = bladeAzimuthSet_rad(:);
bladePitchSet_rad = bladePitchSet_rad(:);

%% -------------------- 2. 机舱质心惯量及塔顶质量矩阵 --------------------
nacelleOffset = [sideOffset;frontOffset;heightOffset];
% Yaw已绕塔筒偏航轴给出，先减去质心侧向/前后偏置产生的移轴项。
% 随后与Noding、Rolling一起统一移轴至塔顶，避免重复计入Yaw移轴项。
yawInertia_aboutMassCenter = yawInertia_aboutTowerAxis ...
    - nacelleMass*(sideOffset^2+frontOffset^2);
nacelleInertiaAtMassCenter = diag([nodingInertia_aboutMassCenter, ...
    rollingInertia_aboutMassCenter,yawInertia_aboutMassCenter]);
[nacelleMatrix,nacelleInertiaAtTowerTop] = rigid_mass_at_tower_top( ...
    nacelleMass,nacelleInertiaAtMassCenter,nacelleOffset);

nacelleOut = struct();
nacelleOut.mass = nacelleMass;
nacelleOut.massCenter_xyz = towerTop_xyz+nacelleOffset;
nacelleOut.rTowerTopToMassCenter = nacelleOffset;
nacelleOut.yawInertia_aboutMassCenter = yawInertia_aboutMassCenter;
nacelleOut.inertiaAtMassCenter_global = nacelleInertiaAtMassCenter;
nacelleOut.inertiaAtTowerTop_global = nacelleInertiaAtTowerTop;
nacelleOut.massMatrix6 = nacelleMatrix;

%% -------------------- 3. 风轮参考点及主轴方向 --------------------
switch rotorPosition
    case "upwind"
        rotorPositionSign = 1;
    case "downwind"
        rotorPositionSign = -1;
    otherwise
        error('rotor_position 必须为 Upwind 或 Downwind。');
end
% H为三支叶片变桨轴交点；overhang本身是水平距离，不再乘cos(tilt)。
rTowerTopToHubReference = [ ...
    lateralOffset;rotorPositionSign*overhang;hubVerticalOffset];
hubReference_xyz = towerTop_xyz+rTowerTopToHubReference;
shaftDirection_global = [ ...
    0;rotorPositionSign*cos(tiltAngle_rad);sin(tiltAngle_rad)];

%% -------------------- 4. 轮毂质心惯量及塔顶质量矩阵 --------------------
rHubReferenceToMassCenter = massCenterOffsetAlongShaft*shaftDirection_global;
hubOffset = rTowerTopToHubReference+rHubReferenceToMassCenter;
% 轴对称惯量：J_G = I_perp*I + (I_shaft-I_perp)*e_s*e_s'。
hubInertiaAtMassCenter = inertiaPerpendicularShaft*eye(3) ...
    + (inertiaShaft-inertiaPerpendicularShaft) ...
    * (shaftDirection_global*shaftDirection_global.');
[hubMatrix,hubInertiaAtTowerTop] = rigid_mass_at_tower_top( ...
    hubMass,hubInertiaAtMassCenter,hubOffset);

hubOut = struct();
hubOut.mass = hubMass;
hubOut.rotorPosition = rotorPosition;
hubOut.tiltAngle_rad = tiltAngle_rad;
hubOut.tiltAngle_deg = rad2deg(tiltAngle_rad);
hubOut.coneAngle_rad = coneAngle_rad;
hubOut.shaftDirection_global = shaftDirection_global;
hubOut.hubReference_xyz = hubReference_xyz;
hubOut.massCenter_xyz = towerTop_xyz+hubOffset;
hubOut.rTowerTopToHubReference = rTowerTopToHubReference;
hubOut.rHubReferenceToMassCenter = rHubReferenceToMassCenter;
hubOut.rTowerTopToMassCenter = hubOffset;
hubOut.inertiaAtMassCenter_global = hubInertiaAtMassCenter;
hubOut.inertiaAtTowerTop_global = hubInertiaAtTowerTop;
hubOut.massMatrix6 = hubMatrix;

%% -------------------- 5. 由质量矩求单叶片质心和惯量 --------------------
% Q1=int(s dm)，Q2=int(s^2 dm)，无需再次沿叶片积分。
% 直线质量分布：sG=Q1/m，两个质心横向惯量为Q2-Q1^2/m，展向惯量为0。
massCenterDistanceFromRoot = bladeFirstMassMoment/bladeMass;
inertiaTransverse_atMassCenter = bladeSecondMassMoment ...
    - bladeFirstMassMoment^2/bladeMass;
bladeInertiaAtMassCenter_root = diag([ ...
    inertiaTransverse_atMassCenter,inertiaTransverse_atMassCenter,0]);
bladeMassCenter_root_xyz = [0;0;massCenterDistanceFromRoot];

%% -------------------- 6. 叶根基准坐标系及锥角旋转 --------------------
xRoot_global = -shaftDirection_global;
zRoot0_global = [0;0;1]-xRoot_global*(xRoot_global.'*[0;0;1]);
zRoot0_global = zRoot0_global/norm(zRoot0_global);
yRoot0_global = cross(zRoot0_global,xRoot_global);
rotZeroAzimuth = [xRoot_global,yRoot0_global,zRoot0_global];
% Rotor正锥角向下风。Upwind时本地正锥变朝塔筒，Downwind时反号。
localConeAngle = rotorPositionSign*coneAngle_rad;
coneCos = cos(localConeAngle);
coneSin = sin(localConeAngle);
rotCone = [coneCos,0,coneSin;0,1,0;-coneSin,0,coneCos];

%% -------------------- 7. 三支叶片逐支等效到塔顶 --------------------
massMatrix6_byBlade = zeros(6,6,3);
massCenter_xyz_byBlade = zeros(3,3);
rTowerTopToMassCenter_byBlade = zeros(3,3);
bladeResult = struct([]);
for bladeIdx = 1:3
    azimuthCos = cos(bladeAzimuthSet_rad(bladeIdx));
    azimuthSin = sin(bladeAzimuthSet_rad(bladeIdx));
    rotAzimuth = [1,0,0;0,azimuthCos,-azimuthSin;0,azimuthSin,azimuthCos];
    pitchCos = cos(bladePitchSet_rad(bladeIdx));
    pitchSin = sin(bladePitchSet_rad(bladeIdx));
    rotPitch = [pitchCos,pitchSin,0;-pitchSin,pitchCos,0;0,0,1];
    rotRootToGlobal = rotZeroAzimuth*rotAzimuth*rotCone*rotPitch;

    % A -> H -> 叶根 -> 叶片质心。L只改变位置，不改变质量。
    rTowerTopToBladeRoot = rTowerTopToHubReference ...
        + rotRootToGlobal*[0;0;bladeRootLength];
    bladeOffset = rTowerTopToBladeRoot ...
        + rotRootToGlobal*bladeMassCenter_root_xyz;
    bladeInertiaAtMassCenter = rotRootToGlobal ...
        * bladeInertiaAtMassCenter_root*rotRootToGlobal.';
    [bladeMatrix,bladeInertiaAtTowerTop] = rigid_mass_at_tower_top( ...
        bladeMass,bladeInertiaAtMassCenter,bladeOffset);
    massMatrix6_byBlade(:,:,bladeIdx) = bladeMatrix;
    rTowerTopToMassCenter_byBlade(:,bladeIdx) = bladeOffset;
    massCenter_xyz_byBlade(:,bladeIdx) = towerTop_xyz+bladeOffset;

    % 关于通过H的风轮主轴惯量，供与Bladed单叶片汇总值比较。
    skewHubOffset = skew_matrix(bladeOffset-rTowerTopToHubReference);
    bladeInertiaAtHubReference = bladeInertiaAtMassCenter ...
        + bladeMass*(skewHubOffset.'*skewHubOffset);
    bladeResult(bladeIdx).mass = bladeMass;
    bladeResult(bladeIdx).azimuth_rad = bladeAzimuthSet_rad(bladeIdx);
    bladeResult(bladeIdx).pitch_rad = bladePitchSet_rad(bladeIdx);
    bladeResult(bladeIdx).rotRootToGlobal = rotRootToGlobal;
    bladeResult(bladeIdx).rootReference_xyz = towerTop_xyz+rTowerTopToBladeRoot;
    bladeResult(bladeIdx).massCenter_xyz = massCenter_xyz_byBlade(:,bladeIdx);
    bladeResult(bladeIdx).rTowerTopToMassCenter = bladeOffset;
    bladeResult(bladeIdx).inertiaAtMassCenter_global = bladeInertiaAtMassCenter;
    bladeResult(bladeIdx).inertiaAtTowerTop_global = bladeInertiaAtTowerTop;
    bladeResult(bladeIdx).inertiaAboutShaft = shaftDirection_global.' ...
        * bladeInertiaAtHubReference*shaftDirection_global;
    bladeResult(bladeIdx).massMatrix6 = bladeMatrix;
end

threeBladesOut = struct();
threeBladesOut.mass = 3*bladeMass;
threeBladesOut.massCenter_xyz = towerTop_xyz+mean(rTowerTopToMassCenter_byBlade,2);
threeBladesOut.massMatrix6 = sum(massMatrix6_byBlade,3);
threeBladesOut.massMatrix6_byBlade = massMatrix6_byBlade;
threeBladesOut.massCenter_xyz_byBlade = massCenter_xyz_byBlade;
threeBladesOut.rTowerTopToMassCenter = mean(rTowerTopToMassCenter_byBlade,2);
threeBladesOut.inertiaAtTowerTop_global = threeBladesOut.massMatrix6(4:6,4:6);
threeBladesOut.hubReference_xyz = hubReference_xyz;
threeBladesOut.shaftDirection_global = shaftDirection_global;
threeBladesOut.bladeRootLength = bladeRootLength;
threeBladesOut.coneAngle_rad = coneAngle_rad;
threeBladesOut.tiltAngle_rad = tiltAngle_rad;
threeBladesOut.singleBlade.mass = bladeMass;
threeBladesOut.singleBlade.firstMassMoment = bladeFirstMassMoment;
threeBladesOut.singleBlade.secondMassMoment = bladeSecondMassMoment;
threeBladesOut.singleBlade.massCenter_root_xyz = bladeMassCenter_root_xyz;
threeBladesOut.singleBlade.inertiaAtMassCenter_root = bladeInertiaAtMassCenter_root;
threeBladesOut.blade = bladeResult;

%% -------------------- 8. 完整RNA矩阵、总质心及总惯量 --------------------
% 各部件共享塔顶的平动和转动速度；总动能可加，故完整质量矩阵可加。
massMatrix6 = nacelleMatrix+hubMatrix+threeBladesOut.massMatrix6;
rnaMass = nacelleMass+hubMass+3*bladeMass;
% 使用相对塔顶的偏置求总质心，避免绝对坐标的相减误差。
rTowerTopToMassCenter = (nacelleMass*nacelleOffset+hubMass*hubOffset ...
    + bladeMass*sum(rTowerTopToMassCenter_byBlade,2))/rnaMass;
massCenter_xyz = towerTop_xyz+rTowerTopToMassCenter;
inertiaAtTowerTop_global = massMatrix6(4:6,4:6);
skewTotalOffset = skew_matrix(rTowerTopToMassCenter);
inertiaAtMassCenter_global = inertiaAtTowerTop_global ...
    - rnaMass*(skewTotalOffset.'*skewTotalOffset);

rnaOut = struct();
rnaOut.towerTop_xyz = towerTop_xyz;
rnaOut.mass = rnaMass;
rnaOut.massCenter_xyz = massCenter_xyz;
rnaOut.rTowerTopToMassCenter = rTowerTopToMassCenter;
rnaOut.inertiaAtTowerTop_global = inertiaAtTowerTop_global;
rnaOut.inertiaAtMassCenter_global = inertiaAtMassCenter_global;
rnaOut.massMatrix6 = massMatrix6;
rnaOut.nacelle = nacelleOut;
rnaOut.hub = hubOut;
rnaOut.threeBlades = threeBladesOut;
end

function [massMatrix6,inertiaAtTowerTop_global] = rigid_mass_at_tower_top( ...
    mass,inertiaAtMassCenter_global,rTowerTopToMassCenter)
% vG=vA-S(rAG)*omega。由刚体动能构造塔顶平动、转动及耦合分块。
skewOffset = skew_matrix(rTowerTopToMassCenter);
inertiaAtTowerTop_global = inertiaAtMassCenter_global ...
    + mass*(skewOffset.'*skewOffset);
massMatrix6 = [mass*eye(3),-mass*skewOffset; ...
    mass*skewOffset,inertiaAtTowerTop_global];
end

function S = skew_matrix(r)
% S(r)*a = r×a；用于质心速度映射及惯量平行移轴。
S = [0,-r(3),r(2);r(3),0,-r(1);-r(2),r(1),0];
end
