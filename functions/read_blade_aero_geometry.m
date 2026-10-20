function bladeAeroData = read_blade_aero_geometry(bladeGeoFile)
%READ_BLADE_AERO_GEOMETRY 仅读取气动截面，不构建柔性叶片质量/刚度。
%
% 输入 bladeGeoFile：沿用原 blade_geo.xlsx 的纯数值表，每列一个截面。
%   第1行：Distance along blade，物理展向长度坐标，m，用于积分ds。
%   第2行：Distance along blade root Z-axis，叶根Z坐标，m，用于力臂。
%   第3行：弦长，m；第4行：气动扭角，deg；第5行：厚弦比，%。
%   第10行：Foil section 编号。
% 第1、2行均沿用表中叶根长度基准，不再增加轮毂半径。
%
% 输出各截面量均为列向量。本阶段沿用直线RNA假设，气动力作用点
% 取展向中心线[0;0;z]，不读取第6~9行的预弯、扫掠和截面中心偏置。
% 第1行和第2行用途不同，不用欧氏力臂长度替代物理积分长度。

geometry = readmatrix(bladeGeoFile);
bladeAeroData.distanceAlongBlade = geometry(1, :).';
bladeAeroData.spanPositionFromRoot = geometry(2, :).';
bladeAeroData.chord = geometry(3, :).';
bladeAeroData.aeroTwist_rad = deg2rad(geometry(4, :).');
bladeAeroData.thickness_pct = geometry(5, :).';
bladeAeroData.foilSection = geometry(10, :).';

% 不排序，以保留Excel截面与翼型信息的对应关系。
if numel(bladeAeroData.distanceAlongBlade) < 2 || ...
        any(diff(bladeAeroData.distanceAlongBlade) <= 0)
    error('blade_geo.xlsx 第1行需要至少两个严格递增的展向位置。');
end
end
