function aerofoilSectionSet = build_aerofoil_section_set(aerofoilDatabase, sectionFile)
%BUILD_AEROFOIL_SECTION_SET
% 根据外部文件构造翼型截面插值规则。
%
% 输入：
%   aerofoilDatabase
%       build_aerofoil_database 生成的翼型极曲线数据库。
%
%   sectionFile
%       外部翼型 Section 规则文件，例如：
%       'AerofoilSectionSet.xlsx'
%
%       文件中应包含 Sheet: Sections
%       并包含两列：
%           SectionNumber
%           FoilName
%
% 输出：
%   aerofoilSectionSet
%       翼型截面插值规则结构体数组。
%
% 外部文件示例：
%   SectionNumber    FoilName
%   1                MY-B-40_VG01
%   1                MY-B-100
%   2                MY-B-35
%   2                MY-B-40
%   3                MY-B-30
%   3                MY-B-35

%% -------------------- 读取外部 Section 文件 --------------------
sectionTable = readtable(sectionFile, ...
    'Sheet', 'Sections', ...
    'VariableNamingRule', 'preserve');

sectionNumberData = sectionTable.SectionNumber;
foilNameData      = string(sectionTable.FoilName);

sectionNumberList = unique(sectionNumberData, 'stable');
numSection = numel(sectionNumberList);

aerofoilSectionSet = struct([]);

%% -------------------- 构造 Section Set --------------------
for sectionIdx = 1:numSection

    sectionNumber = sectionNumberList(sectionIdx);

    sectionRowIdx = sectionNumberData == sectionNumber;

    foilNameList = foilNameData(sectionRowIdx);
    foilNameList = foilNameList(:);

    numFoil = numel(foilNameList);
    tcPercentList = zeros(numFoil, 1);

    for foilIdx = 1:numFoil
        tcPercentList(foilIdx) = get_aerofoil_tc_percent( ...
            aerofoilDatabase, foilNameList(foilIdx));
    end

    % 按厚度从小到大排序，保证后续插值方向明确
    [tcPercentList, sortIdx] = sort(tcPercentList);
    foilNameList = foilNameList(sortIdx);

    if any(diff(tcPercentList) <= 0)
        error(['SectionNumber %g 中存在重复或无效的翼型厚弦比，' ...
            '无法按厚弦比插值。'], sectionNumber);
    end

    aerofoilSectionSet(sectionIdx).sectionNumber = sectionNumber;

    aerofoilSectionSet(sectionIdx).foilNameList = foilNameList;
    aerofoilSectionSet(sectionIdx).tcPercentList = tcPercentList;

    % 下面几个字段保留，是为了兼容你之前的 query_aerofoil_section 写法
    aerofoilSectionSet(sectionIdx).lowerFoilName = foilNameList(1);
    aerofoilSectionSet(sectionIdx).upperFoilName = foilNameList(end);

    aerofoilSectionSet(sectionIdx).lowerTcPercent = tcPercentList(1);
    aerofoilSectionSet(sectionIdx).upperTcPercent = tcPercentList(end);
end

end


%% ========================================================================
function tc_percent = get_aerofoil_tc_percent(aerofoilDatabase, aerofoilName)
%GET_AEROFOIL_TC_PERCENT
% 根据翼型名称从 aerofoilDatabase 中提取厚弦比百分数。

nameList = string({aerofoilDatabase.foil.name});
foilIdx = find(nameList == string(aerofoilName), 1);

if isempty(foilIdx)
    error('翼型数据库中未找到 %s。', char(string(aerofoilName)));
end

tc_percent = aerofoilDatabase.foil(foilIdx).tc_percent;

end
