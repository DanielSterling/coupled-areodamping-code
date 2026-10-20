function [NodeData_FE, ElemData, Conn_FE, Tatong_sorted, StationData] = ...
    tower_elem_props_from_Tatong_gauss_FE(Tatong)
%==========================================================================
% 功能：
% 1) 按高度对塔筒 station 排序
% 2) 必要时补充 z = 0 station
% 3) 处理同高 station：
%    - 同高且属性相同：去重
%    - 同高但属性不同：保留，用于描述截面突变
% 4) 同一高度只建立一个有限元节点
% 5) 在相邻不同高度之间建立空间梁单元
% 6) 保存单元两端几何属性和力学属性
%
% 输入：
%   Tatong : 每列对应一个 station
%
% 数据行：
%   row 1 : Height above surface                  [m]
%   row 2 : Outside diameter                      [m]
%   row 3 : Mass per unit length                  [kg/m]
%   row 4 : Bending stiffness                     [N*m^2]
%   row 5 : Shear stiffness                       [N]
%   row 6 : Torsional stiffness                   [N*m^2]
%   row 7 : Mass moment of inertia per unit length
%   row 8 : Wall thickness                        [mm]
%
% 输出：
%   NodeData_FE   : 有限元节点信息
%   ElemData      : 单元信息
%   Conn_FE       : 单元连接矩阵
%   Tatong_sorted : 排序、补零和去重后的 station 数据
%   StationData   : station 层面的几何与力学属性
%
% 同高不同属性 station 的顺序约定：
%   第一个 station 供该高度以下的单元使用；
%   第二个 station 供该高度以上的单元使用。
%
% 本函数不执行：
%   - 高斯积分
%   - 单元刚度矩阵计算
%   - 单元质量矩阵计算
%   - 边界条件处理
%==========================================================================

    %% -------------------- 输入数据行定义 --------------------

    row_z  = 1;
    row_D  = 2;
    row_m  = 3;
    row_EI = 4;
    row_GA = 5;
    row_GJ = 6;
    row_Ip = 7;
    row_t  = 8;

    tol_z    = 1e-10;
    tol_prop = 1e-10;

    propertyRows = ...
        [row_D, row_m, row_EI, row_GA, row_GJ, row_Ip, row_t];

    if isempty(Tatong) || size(Tatong, 1) < 8
        error('Tatong 至少需要包含 8 行塔筒 station 数据。');
    end

    if any(~isfinite(Tatong(:)))
        error('Tatong 中存在 NaN 或 Inf。');
    end


    %% -------------------- 1. 按高度排序 --------------------

    [~, sortId] = sort(Tatong(row_z, :), 'ascend');

    Tatong_sorted = Tatong(:, sortId);


    %% -------------------- 2. 补充 z = 0 station --------------------

    baseStationAdded = false;

    z_min = Tatong_sorted(row_z, 1);

    if z_min < -tol_z

        error( ...
            '检测到负的塔筒 station 高度：z_min = %.6g m。', ...
            z_min);

    elseif z_min > tol_z

        newStation = Tatong_sorted(:, 1);
        newStation(row_z) = 0;

        Tatong_sorted = [newStation, Tatong_sorted];

        baseStationAdded = true;

    else

        % 将数值上非常接近 0 的高度统一设置为 0
        nearZero = abs(Tatong_sorted(row_z, :)) <= tol_z;
        Tatong_sorted(row_z, nearZero) = 0;

    end


    %% -------------------- 3. 同高 station 去重 --------------------

    Tatong_clean = zeros(size(Tatong_sorted, 1), 0);

    nInputStation = size(Tatong_sorted, 2);

    i = 1;

    while i <= nInputStation

        zGroup = Tatong_sorted(row_z, i);

        % 找出当前同高 station 组
        j = i;

        while j <= nInputStation && ...
                abs(Tatong_sorted(row_z, j) - zGroup) <= tol_z

            j = j + 1;

        end

        groupId = i:(j - 1);

        % 保存当前高度下互不重复的 station
        retainedId = zeros(1, 0);

        for k = groupId

            currentProperty = Tatong_sorted(propertyRows, k);

            isDuplicate = false;

            for r = retainedId

                retainedProperty = Tatong_sorted(propertyRows, r);

                propertyScale = max( ...
                    [ones(numel(propertyRows), 1), ...
                     abs(currentProperty), ...
                     abs(retainedProperty)], ...
                    [], ...
                    2);

                propertyDifference = ...
                    abs(currentProperty - retainedProperty);

                if all(propertyDifference <= ...
                        tol_prop .* propertyScale)

                    isDuplicate = true;
                    break;

                end
            end

            if ~isDuplicate
                retainedId(end + 1) = k; %#ok<AGROW>
            end
        end

        % 一个高度最多允许上下两组截面属性
        if numel(retainedId) > 2

            error( ...
                ['高度 z = %.6g m 处存在超过两组不同的截面属性，' ...
                 '无法判断各组属性分别属于哪个单元。'], ...
                zGroup);
        end

        groupData = Tatong_sorted(:, retainedId);

        % 将容差范围内的高度统一成同一个数值
        groupData(row_z, :) = zGroup;

        Tatong_clean = [Tatong_clean, groupData]; %#ok<AGROW>

        i = j;
    end

    Tatong_sorted = Tatong_clean;


    %% -------------------- 4. 构造 StationData --------------------

    nSta = size(Tatong_sorted, 2);

    StationData.id = (1:nSta)';

    StationData.z = Tatong_sorted(row_z, :).';
    StationData.x = zeros(nSta, 1);
    StationData.y = zeros(nSta, 1);

    StationData.xyz = ...
        [StationData.x, StationData.y, StationData.z];

    % 原始输入属性
    StationData.D = Tatong_sorted(row_D, :).';
    StationData.t = Tatong_sorted(row_t, :).' / 1000;

    StationData.m = Tatong_sorted(row_m, :).';

    StationData.EIy = Tatong_sorted(row_EI, :).';
    StationData.EIz = Tatong_sorted(row_EI, :).';

    StationData.GAy = Tatong_sorted(row_GA, :).';
    StationData.GAz = Tatong_sorted(row_GA, :).';

    StationData.GJ = Tatong_sorted(row_GJ, :).';

    % 单位长度质量转动惯量
    StationData.Ip = Tatong_sorted(row_Ip, :).';

    StationData.baseStationAdded = baseStationAdded;


    %% -------------------- 5. 圆管截面几何属性 --------------------

    StationData.Di = StationData.D - 2 .* StationData.t;

    if any(StationData.D <= 0) || ...
            any(StationData.t <= 0) || ...
            any(StationData.Di <= 0)

        error('塔筒外径或壁厚数据不合理。');
    end

    if any(StationData.m <= 0)
        error('塔筒单位长度质量必须大于 0。');
    end

    % 截面积
    StationData.A = pi / 4 .* ...
        (StationData.D.^2 - StationData.Di.^2);

    % 截面惯性矩
    StationData.I = pi / 64 .* ...
        (StationData.D.^4 - StationData.Di.^4);

    % 几何极惯性矩
    StationData.J = pi / 32 .* ...
        (StationData.D.^4 - StationData.Di.^4);

    % 根据 EI 和几何惯性矩反算等效弹性模量
    StationData.E = StationData.EIy ./ StationData.I;

    % 塔筒轴向刚度
    StationData.EA = StationData.E .* StationData.A;


    %% -------------------- 6. 建立唯一 FE 节点 --------------------

    z_sta = StationData.z;

    ic = zeros(nSta, 1);
    z_unique = zeros(nSta, 1);

    nNode = 1;

    z_unique(1) = z_sta(1);
    ic(1) = 1;

    for i = 2:nSta

        if abs(z_sta(i) - z_unique(nNode)) <= tol_z

            ic(i) = nNode;

        else

            nNode = nNode + 1;

            z_unique(nNode) = z_sta(i);
            ic(i) = nNode;

        end
    end

    z_unique = z_unique(1:nNode);

    NodeData_FE.id = (1:nNode)';
    NodeData_FE.x = zeros(nNode, 1);
    NodeData_FE.y = zeros(nNode, 1);
    NodeData_FE.z = z_unique;

    NodeData_FE.xyz = ...
        [NodeData_FE.x, NodeData_FE.y, NodeData_FE.z];


    %% -------------------- 7. 建立有限元单元 --------------------

    ElemData = struct([]);

    Conn_FE = zeros(max(nSta - 1, 0), 2);

    nElem = 0;

    for i = 1:(nSta - 1)

        % 同高 station 之间不建立单元
        if abs(StationData.z(i + 1) - StationData.z(i)) <= tol_z
            continue;
        end

        nElem = nElem + 1;

        station1 = i;
        station2 = i + 1;

        node1 = ic(station1);
        node2 = ic(station2);

        xyz1 = NodeData_FE.xyz(node1, :);
        xyz2 = NodeData_FE.xyz(node2, :);

        elementVector = xyz2 - xyz1;
        elementLength = norm(elementVector);

        Conn_FE(nElem, :) = [node1, node2];


        %% -------------------- 拓扑信息 --------------------

        ElemData(nElem).id = nElem;

        ElemData(nElem).node1 = node1;
        ElemData(nElem).node2 = node2;
        ElemData(nElem).nodes = [node1, node2];

        ElemData(nElem).station1 = station1;
        ElemData(nElem).station2 = station2;
        ElemData(nElem).stations = [station1, station2];


        %% -------------------- 几何信息 --------------------

        ElemData(nElem).xyz1 = xyz1;
        ElemData(nElem).xyz2 = xyz2;

        ElemData(nElem).x1 = xyz1(1);
        ElemData(nElem).y1 = xyz1(2);
        ElemData(nElem).z1 = xyz1(3);

        ElemData(nElem).x2 = xyz2(1);
        ElemData(nElem).y2 = xyz2(2);
        ElemData(nElem).z2 = xyz2(3);

        ElemData(nElem).Le = elementLength;

        % 单元轴线在全局坐标系中的单位方向
        ElemData(nElem).axisVector = ...
            elementVector / elementLength;


        %% -------------------- 截面尺寸 --------------------

        ElemData(nElem).D1 = StationData.D(station1);
        ElemData(nElem).D2 = StationData.D(station2);

        ElemData(nElem).t1 = StationData.t(station1);
        ElemData(nElem).t2 = StationData.t(station2);

        ElemData(nElem).Dmean = ...
            0.5 * (ElemData(nElem).D1 + ElemData(nElem).D2);


        %% -------------------- 质量属性 --------------------

        ElemData(nElem).m1 = StationData.m(station1);
        ElemData(nElem).m2 = StationData.m(station2);

        ElemData(nElem).Ip1 = StationData.Ip(station1);
        ElemData(nElem).Ip2 = StationData.Ip(station2);


        %% -------------------- 刚度属性 --------------------

        ElemData(nElem).EA1 = StationData.EA(station1);
        ElemData(nElem).EA2 = StationData.EA(station2);

        ElemData(nElem).EIy1 = StationData.EIy(station1);
        ElemData(nElem).EIy2 = StationData.EIy(station2);

        ElemData(nElem).EIz1 = StationData.EIz(station1);
        ElemData(nElem).EIz2 = StationData.EIz(station2);

        ElemData(nElem).GAy1 = StationData.GAy(station1);
        ElemData(nElem).GAy2 = StationData.GAy(station2);

        ElemData(nElem).GAz1 = StationData.GAz(station1);
        ElemData(nElem).GAz2 = StationData.GAz(station2);

        ElemData(nElem).GJ1 = StationData.GJ(station1);
        ElemData(nElem).GJ2 = StationData.GJ(station2);


        %% -------------------- 元信息 --------------------

        ElemData(nElem).sectionType = 'tube';

    end


    %% -------------------- 8. 裁剪连接矩阵 --------------------

    Conn_FE = Conn_FE(1:nElem, :);

end