function [Cl, Cd, Cm, dCl_dalpha, dCd_dalpha, dCm_dalpha] = ...
    query_aerofoil_section( ...
    aerofoilDatabase, aerofoilSectionSet, ...
    foilSectionNumber, tc_percent, alpha_rad)
%QUERY_AEROFOIL_SECTION
% 根据 foil section 编号、厚弦比和攻角查询截面气动系数。
%
% 输入：
%   aerofoilDatabase
%       build_aerofoil_database 生成的翼型极曲线数据库
%
%   aerofoilSectionSet
%       build_aerofoil_section_set 生成的翼型截面插值规则
%
%   foilSectionNumber
%       foil section 编号，可以是标量，也可以是 nSection×1
%
%   tc_percent
%       截面厚弦比百分数，可以是标量，也可以是 nSection×1
%
%   alpha_rad
%       截面攻角，单位 rad，可以是标量，也可以是 nSection×1
%
% 输出：
%   Cl, Cd, Cm
%       升力系数、阻力系数、俯仰力矩系数
%
%   dCl_dalpha, dCd_dalpha, dCm_dalpha
%       气动系数对攻角的导数，单位 1/rad

%% -------------------- 输入整理 --------------------
numQuery = max([ ...
    numel(foilSectionNumber), ...
    numel(tc_percent), ...
    numel(alpha_rad)]);

inputLengthSet = [ ...
    numel(foilSectionNumber), ...
    numel(tc_percent), ...
    numel(alpha_rad)];

if any(inputLengthSet ~= 1 & inputLengthSet ~= numQuery)
    error(['foilSectionNumber、tc_percent 和 alpha_rad 的长度必须相同，', ...
        '或者其中某项为标量。']);
end

if numel(alpha_rad) > 1
    outputSize = size(alpha_rad);
elseif numel(tc_percent) > 1
    outputSize = size(tc_percent);
elseif numel(foilSectionNumber) > 1
    outputSize = size(foilSectionNumber);
else
    outputSize = size(alpha_rad);
end

foilSectionNumber = foilSectionNumber(:);
tc_percent = tc_percent(:);
alpha_rad = alpha_rad(:);

if numel(foilSectionNumber) == 1
    foilSectionNumber = repmat(foilSectionNumber, numQuery, 1);
end

if numel(tc_percent) == 1
    tc_percent = repmat(tc_percent, numQuery, 1);
end

if numel(alpha_rad) == 1
    alpha_rad = repmat(alpha_rad, numQuery, 1);
end

%% -------------------- 初始化输出 --------------------
Cl = zeros(numQuery, 1);
Cd = zeros(numQuery, 1);
Cm = zeros(numQuery, 1);

dCl_dalpha = zeros(numQuery, 1);
dCd_dalpha = zeros(numQuery, 1);
dCm_dalpha = zeros(numQuery, 1);

sectionNumberList = [aerofoilSectionSet.sectionNumber];

%% -------------------- 逐截面查询 --------------------
for queryIdx = 1:numQuery

    currentSectionNumber = foilSectionNumber(queryIdx);
    currentTcPercent = tc_percent(queryIdx);
    currentAlpha = alpha_rad(queryIdx);

    aerofoilSectionIdx = find(sectionNumberList == currentSectionNumber, 1);

    if isempty(aerofoilSectionIdx)
        error('未找到 Foil section %g 的翼型插值规则。', ...
            currentSectionNumber);
    end

    aerofoilSection = aerofoilSectionSet(aerofoilSectionIdx);

    foilNameList = aerofoilSection.foilNameList;
    tcPercentList = aerofoilSection.tcPercentList;

    % 使用包围当前厚弦比的相邻翼型进行分段线性插值。超出数据库范围
    % 时取最近端点；只有一个翼型时直接查询，避免除以零。
    if numel(foilNameList) == 1 || currentTcPercent <= tcPercentList(1)
        lowerIdx = 1;
        upperIdx = 1;
        weightUpper = 0;
    elseif currentTcPercent >= tcPercentList(end)
        lowerIdx = numel(foilNameList);
        upperIdx = lowerIdx;
        weightUpper = 0;
    else
        upperIdx = find(tcPercentList >= currentTcPercent, 1);
        lowerIdx = upperIdx - 1;
        tcLower = tcPercentList(lowerIdx);
        tcUpper = tcPercentList(upperIdx);
        weightUpper = ...
            (currentTcPercent - tcLower) / (tcUpper - tcLower);
    end

    lowerFoilName = foilNameList(lowerIdx);
    upperFoilName = foilNameList(upperIdx);

    [ClLower, CdLower, CmLower, dClLower, dCdLower, dCmLower] = ...
        aerofoilDatabase.query(lowerFoilName, currentAlpha);

    if upperIdx == lowerIdx
        ClUpper = ClLower;
        CdUpper = CdLower;
        CmUpper = CmLower;
        dClUpper = dClLower;
        dCdUpper = dCdLower;
        dCmUpper = dCmLower;
    else
        [ClUpper, CdUpper, CmUpper, dClUpper, dCdUpper, dCmUpper] = ...
            aerofoilDatabase.query(upperFoilName, currentAlpha);
    end

    Cl(queryIdx) = ...
        (1 - weightUpper) * ClLower + weightUpper * ClUpper;

    Cd(queryIdx) = ...
        (1 - weightUpper) * CdLower + weightUpper * CdUpper;

    Cm(queryIdx) = ...
        (1 - weightUpper) * CmLower + weightUpper * CmUpper;

    dCl_dalpha(queryIdx) = ...
        (1 - weightUpper) * dClLower + weightUpper * dClUpper;

    dCd_dalpha(queryIdx) = ...
        (1 - weightUpper) * dCdLower + weightUpper * dCdUpper;

    dCm_dalpha(queryIdx) = ...
        (1 - weightUpper) * dCmLower + weightUpper * dCmUpper;
end

%% -------------------- 恢复输出尺寸 --------------------
Cl = reshape(Cl, outputSize);
Cd = reshape(Cd, outputSize);
Cm = reshape(Cm, outputSize);

dCl_dalpha = reshape(dCl_dalpha, outputSize);
dCd_dalpha = reshape(dCd_dalpha, outputSize);
dCm_dalpha = reshape(dCm_dalpha, outputSize);

end
