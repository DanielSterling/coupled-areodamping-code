function aerofoilDatabase = build_aerofoil_database(xfile)
%BUILD_AEROFOIL_DATABASE
% 构造多翼型气动极曲线数据库。
%
% 输入：
%   xfile
%       翼型气动数据库 Excel 文件，例如 "AerofoilDatabase.xlsx"
%
% 输出：
%   aerofoilDatabase
%       翼型极曲线数据库结构体
%
% 查询方式：
%   [Cl, Cd, Cm, dCl_dalpha, dCd_dalpha, dCm_dalpha] = ...
%       aerofoilDatabase.query(aerofoilName, alpha_rad)

%% 1. 读取所有 sheet 名称
sheetNameList = string(sheetnames(xfile));
numAerofoil = numel(sheetNameList);

%% 2. 逐个读取翼型极曲线
foil = struct([]);

for foilIdx = 1:numAerofoil

    aerofoilName = sheetNameList(foilIdx);
    sheetName    = sheetNameList(foilIdx);

    data = readmatrix(xfile, "Sheet", sheetName);

    alpha_deg = data(:, 1);
    Cl        = data(:, 2);
    Cd        = data(:, 3);
    Cm        = data(:, 4);

    % 去掉表头或空白行
    validRow = ~isnan(alpha_deg) & ~isnan(Cl) & ~isnan(Cd) & ~isnan(Cm);

    alpha_deg = alpha_deg(validRow);
    Cl        = Cl(validRow);
    Cd        = Cd(validRow);
    Cm        = Cm(validRow);

    % 攻角 degree -> rad
    alpha_rad = deg2rad(alpha_deg);

    % 按攻角排序
    [alpha_rad, sortIdx] = sort(alpha_rad(:));

    alpha_deg = alpha_deg(sortIdx);
    Cl        = Cl(sortIdx);
    Cd        = Cd(sortIdx);
    Cm        = Cm(sortIdx);

    % 插值器：C(alpha)
    FCl = griddedInterpolant(alpha_rad, Cl, "linear", "nearest");
    FCd = griddedInterpolant(alpha_rad, Cd, "linear", "nearest");
    FCm = griddedInterpolant(alpha_rad, Cm, "linear", "nearest");

    % 数值求导步长
    alphaStep = min(diff(alpha_rad)) / 10;

    foil(foilIdx).name       = aerofoilName;
    foil(foilIdx).sheetName  = sheetName;
    foil(foilIdx).tc_percent = parse_tc_percent_from_foil_name(aerofoilName);

    foil(foilIdx).alpha_deg  = alpha_deg;
    foil(foilIdx).alpha_rad  = alpha_rad;

    foil(foilIdx).Cl = Cl;
    foil(foilIdx).Cd = Cd;
    foil(foilIdx).Cm = Cm;

    foil(foilIdx).FCl = FCl;
    foil(foilIdx).FCd = FCd;
    foil(foilIdx).FCm = FCm;

    foil(foilIdx).alphaStep = alphaStep;

    foil(foilIdx).query = @(alphaQuery_rad) query_single_aerofoil( ...
        FCl, FCd, FCm, alphaQuery_rad, alphaStep, ...
        min(alpha_rad), max(alpha_rad));
end

%% 3. 打包数据库
aerofoilDatabase = struct();

aerofoilDatabase.xfile = xfile;
aerofoilDatabase.foil  = foil;

aerofoilDatabase.nameList = string({foil.name}).';
aerofoilDatabase.tcPercentList = [foil.tc_percent].';

aerofoilDatabase.query = @(aerofoilName, alphaQuery_rad) ...
    query_aerofoil_database(foil, aerofoilName, alphaQuery_rad);

end


%% ========================================================================
function [Cl, Cd, Cm, dCl_dalpha, dCd_dalpha, dCm_dalpha] = ...
    query_aerofoil_database(foil, aerofoilName, alpha_rad)
%QUERY_AEROFOIL_DATABASE
% 根据翼型名称查询气动系数和攻角导数。

nameList = string({foil.name});
foilIdx = find(nameList == string(aerofoilName), 1);

[Cl, Cd, Cm, dCl_dalpha, dCd_dalpha, dCm_dalpha] = ...
    foil(foilIdx).query(alpha_rad);

end


%% ========================================================================
function [Cl, Cd, Cm, dCl_dalpha, dCd_dalpha, dCm_dalpha] = ...
    query_single_aerofoil(FCl, FCd, FCm, alpha_rad, alphaStep, alphaMin, alphaMax)
%QUERY_SINGLE_AEROFOIL
% 查询单个翼型的 Cl、Cd、Cm 及其攻角导数。

alphaSize = size(alpha_rad);
alpha_rad = alpha_rad(:);

Cl = FCl(alpha_rad);
Cd = FCd(alpha_rad);
Cm = FCm(alpha_rad);

dCl_dalpha = diff_alpha(FCl, alpha_rad, alphaStep, alphaMin, alphaMax);
dCd_dalpha = diff_alpha(FCd, alpha_rad, alphaStep, alphaMin, alphaMax);
dCm_dalpha = diff_alpha(FCm, alpha_rad, alphaStep, alphaMin, alphaMax);

Cl = reshape(Cl, alphaSize);
Cd = reshape(Cd, alphaSize);
Cm = reshape(Cm, alphaSize);

dCl_dalpha = reshape(dCl_dalpha, alphaSize);
dCd_dalpha = reshape(dCd_dalpha, alphaSize);
dCm_dalpha = reshape(dCm_dalpha, alphaSize);

end


%% ========================================================================
function dC_dalpha = diff_alpha(F, alpha_rad, alphaStep, alphaMin, alphaMax)
%DIFF_ALPHA
% 对插值函数 C(alpha) 沿 alpha 方向做数值求导。

alphaPlus  = alpha_rad + alphaStep;
alphaMinus = alpha_rad - alphaStep;

dC_dalpha = zeros(size(alpha_rad));

% 中心差分
idxMiddle = alphaMinus >= alphaMin & alphaPlus <= alphaMax;

dC_dalpha(idxMiddle) = ...
    (F(alphaPlus(idxMiddle)) - F(alphaMinus(idxMiddle))) / (2 * alphaStep);

% 下边界：前向差分
idxLower = alphaMinus < alphaMin;

dC_dalpha(idxLower) = ...
    (F(alpha_rad(idxLower) + alphaStep) - F(alpha_rad(idxLower))) / alphaStep;

% 上边界：后向差分
idxUpper = alphaPlus > alphaMax;

dC_dalpha(idxUpper) = ...
    (F(alpha_rad(idxUpper)) - F(alpha_rad(idxUpper) - alphaStep)) / alphaStep;

end


%% ========================================================================
function tc_percent = parse_tc_percent_from_foil_name(aerofoilName)
%PARSE_TC_PERCENT_FROM_FOIL_NAME
% 从翼型名称中解析厚度百分数。
%
% 示例：
%   MY-B-40       -> 40
%   MY-B-40_VG01  -> 40

nums = regexp(string(aerofoilName), '\d+\.?\d*', 'match');

tc_percent = str2double(nums{1});

end