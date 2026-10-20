function [Ke_global, Me_global, Ke_local, Me_local, ...
          T_gl, info, R_gl] = ...
    beam3d_tower_TIM_station(elem, opts)
%==========================================================================
% 两节点三维 Timoshenko 塔筒梁单元
%
% 本函数直接使用塔筒 station 给出的等效截面属性：
%   EA、EI、GA、GJ、m、Ip
%
% 不再根据 D、t、E、G、rho 重新计算质量和刚度。
%
% 输入：
%   elem : tower_elem_props_from_Tatong_gauss_FE 生成的单元结构体
%
%          必要字段：
%          elem.xyz1, elem.xyz2
%
%          elem.EA1,  elem.EA2
%          elem.EIy1, elem.EIy2
%          elem.EIz1, elem.EIz2
%          elem.GAy1, elem.GAy2
%          elem.GAz1, elem.GAz2
%          elem.GJ1,  elem.GJ2
%
%          elem.m1,   elem.m2
%          elem.Ip1,  elem.Ip2
%
%   opts : 可选计算设置
%          opts.refVec : 局部坐标系参考方向，默认 [0 1 0]
%          opts.R_gl   : 可直接指定局部到全局旋转矩阵
%          opts.nGP    : 刚度和质量积分点数，默认 2
%
% 局部自由度顺序：
%   [u1 v1 w1 rx1 ry1 rz1 u2 v2 w2 rx2 ry2 rz2]
%
% 局部坐标：
%   local x : 梁轴线方向
%   local y、z : 截面横向
%
% 坐标变换：
%   d_global = T_gl * d_local
%
%   Ke_global = T_gl * Ke_local * T_gl'
%   Me_global = T_gl * Me_local * T_gl'
%
% 数值积分：
%   轴向、弯曲、扭转：标准 Gauss 积分
%   剪切项：一点选择性减缩积分，避免细长塔筒剪切锁死
%   质量矩阵：标准 Gauss 积分
%==========================================================================

    %% -------------------- 默认设置 --------------------

    if nargin < 2 || isempty(opts)
        opts = struct();
    end

    if ~isfield(opts, 'refVec')
        opts.refVec = [0 1 0];
    end

    if ~isfield(opts, 'nGP')
        opts.nGP = 2;
    end


    %% -------------------- 单元几何 --------------------

    xyz1 = elem.xyz1(:).';
    xyz2 = elem.xyz2(:).';

    elementVector = xyz2 - xyz1;
    L = norm(elementVector);

    if L <= 0
        error('梁单元长度必须大于 0。');
    end


    %% -------------------- 局部坐标系 --------------------

    if isfield(opts, 'R_gl') && ~isempty(opts.R_gl)

        R_gl = opts.R_gl;

        ex = R_gl(:,1);
        ey = R_gl(:,2);
        ez = R_gl(:,3);

    else

        [R_gl, ex, ey, ez] = ...
            build_local_frame(xyz1, xyz2, opts.refVec);

    end

    T_gl = build_T_matrix(R_gl);


    %% -------------------- 初始化单元矩阵 --------------------

    Ke_local = zeros(12,12);
    Me_local = zeros(12,12);

    % 两节点线性形函数对物理坐标的导数
    dNdx = [-1/L, 1/L];

    % 自然坐标到物理坐标的 Jacobian
    Jmap = L/2;


    %% ====================================================================
    % 1. 轴向、弯曲和扭转刚度
    %% ====================================================================

    [xiSet, weightSet] = gauss_rule_1d(opts.nGP);

    for gp = 1:numel(xiSet)

        xi = xiSet(gp);
        weight = weightSet(gp);

        N1 = (1 - xi)/2;
        N2 = (1 + xi)/2;

        N = [N1, N2];


        %% -------------------- 高斯点属性 --------------------

        EA_gp = ...
            N1*elem.EA1 + N2*elem.EA2;

        EIy_gp = ...
            N1*elem.EIy1 + N2*elem.EIy2;

        EIz_gp = ...
            N1*elem.EIz1 + N2*elem.EIz2;

        GJ_gp = ...
            N1*elem.GJ1 + N2*elem.GJ2;


        %% -------------------- 轴向应变 --------------------
        % epsilon_x = du/dx

        B_axial = zeros(1,12);
        B_axial([1,7]) = dNdx;

        Ke_local = Ke_local + ...
            B_axial.' * EA_gp * B_axial * Jmap * weight;


        %% -------------------- 扭转应变 --------------------
        % kappa_x = d(rx)/dx

        B_torsion = zeros(1,12);
        B_torsion([4,10]) = dNdx;

        Ke_local = Ke_local + ...
            B_torsion.' * GJ_gp * B_torsion * Jmap * weight;


        %% -------------------- 绕 local y 轴弯曲 --------------------
        % w-ry 子系统
        % kappa_y = d(ry)/dx

        B_bending_y = zeros(1,12);
        B_bending_y([5,11]) = dNdx;

        Ke_local = Ke_local + ...
            B_bending_y.' * EIy_gp * B_bending_y * ...
            Jmap * weight;


        %% -------------------- 绕 local z 轴弯曲 --------------------
        % v-rz 子系统
        % kappa_z = d(rz)/dx

        B_bending_z = zeros(1,12);
        B_bending_z([6,12]) = dNdx;

        Ke_local = Ke_local + ...
            B_bending_z.' * EIz_gp * B_bending_z * ...
            Jmap * weight;

    end


    %% ====================================================================
    % 2. 剪切刚度：一点选择性减缩积分
    %
    % 对于细长塔筒，如果对剪切项采用完全积分，
    % 二节点线性 Timoshenko 梁可能发生剪切锁死。
    %% ====================================================================

    xi = 0;
    weight = 2;

    N1 = 0.5;
    N2 = 0.5;

    N = [N1, N2];

    GAy_gp = ...
        N1*elem.GAy1 + N2*elem.GAy2;

    GAz_gp = ...
        N1*elem.GAz1 + N2*elem.GAz2;


    %% -------------------- local y 方向剪切 --------------------
    % v-rz 子系统
    %
    % gamma_y = dv/dx - rz

    B_shear_y = zeros(1,12);

    B_shear_y([2,8])  = dNdx;
    B_shear_y([6,12]) = -N;

    Ke_local = Ke_local + ...
        B_shear_y.' * GAy_gp * B_shear_y * ...
        Jmap * weight;


    %% -------------------- local z 方向剪切 --------------------
    % w-ry 子系统
    %
    % gamma_z = dw/dx + ry

    B_shear_z = zeros(1,12);

    B_shear_z([3,9])  = dNdx;
    B_shear_z([5,11]) = N;

    Ke_local = Ke_local + ...
        B_shear_z.' * GAz_gp * B_shear_z * ...
        Jmap * weight;


    %% ====================================================================
    % 3. 一致质量矩阵
    %
    % Ip 解释为绕塔筒轴线的单位长度质量极惯量。
    %
    % 对于轴对称圆管：
    %   Ixx_mass = Ip
    %   Iyy_mass = Ip/2
    %   Izz_mass = Ip/2
    %% ====================================================================

    for gp = 1:numel(xiSet)

        xi = xiSet(gp);
        weight = weightSet(gp);

        N1 = (1 - xi)/2;
        N2 = (1 + xi)/2;

        N = [N1, N2];


        %% -------------------- 高斯点质量属性 --------------------

        m_gp = ...
            N1*elem.m1 + N2*elem.m2;

        Ip_gp = ...
            N1*elem.Ip1 + N2*elem.Ip2;

        Ixx_mass_gp = Ip_gp;
        Iyy_mass_gp = 0.5*Ip_gp;
        Izz_mass_gp = 0.5*Ip_gp;


        %% -------------------- 两节点质量子矩阵 --------------------

        NN = N.' * N;

        M_translation = ...
            m_gp * NN * Jmap * weight;

        M_rotation_x = ...
            Ixx_mass_gp * NN * Jmap * weight;

        M_rotation_y = ...
            Iyy_mass_gp * NN * Jmap * weight;

        M_rotation_z = ...
            Izz_mass_gp * NN * Jmap * weight;


        %% -------------------- 三个平动方向 --------------------

        Me_local([1,7],[1,7]) = ...
            Me_local([1,7],[1,7]) + M_translation;

        Me_local([2,8],[2,8]) = ...
            Me_local([2,8],[2,8]) + M_translation;

        Me_local([3,9],[3,9]) = ...
            Me_local([3,9],[3,9]) + M_translation;


        %% -------------------- 三个转动方向 --------------------

        Me_local([4,10],[4,10]) = ...
            Me_local([4,10],[4,10]) + M_rotation_x;

        Me_local([5,11],[5,11]) = ...
            Me_local([5,11],[5,11]) + M_rotation_y;

        Me_local([6,12],[6,12]) = ...
            Me_local([6,12],[6,12]) + M_rotation_z;

    end


    %% -------------------- 数值对称化 --------------------

    Ke_local = 0.5 * (Ke_local + Ke_local.');
    Me_local = 0.5 * (Me_local + Me_local.');


    %% -------------------- 转换到全局坐标系 --------------------

    Ke_global = ...
        T_gl * Ke_local * T_gl.';

    Me_global = ...
        T_gl * Me_local * T_gl.';

    Ke_global = ...
        0.5 * (Ke_global + Ke_global.');

    Me_global = ...
        0.5 * (Me_global + Me_global.');


    %% -------------------- 输出信息 --------------------

    info = struct();

    info.L = L;

    info.ex = ex;
    info.ey = ey;
    info.ez = ez;

    info.R_gl = R_gl;

    info.EA_mean = ...
        0.5*(elem.EA1 + elem.EA2);

    info.EIy_mean = ...
        0.5*(elem.EIy1 + elem.EIy2);

    info.EIz_mean = ...
        0.5*(elem.EIz1 + elem.EIz2);

    info.GAy_mean = ...
        0.5*(elem.GAy1 + elem.GAy2);

    info.GAz_mean = ...
        0.5*(elem.GAz1 + elem.GAz2);

    info.GJ_mean = ...
        0.5*(elem.GJ1 + elem.GJ2);

    info.m_mean = ...
        0.5*(elem.m1 + elem.m2);

    info.Ip_mean = ...
        0.5*(elem.Ip1 + elem.Ip2);

    info.propertySource = ...
        'station equivalent properties';

    info.shearIntegration = ...
        'one-point reduced integration';

end


%% ========================================================================
% 局部函数：建立梁单元局部坐标系
%% ========================================================================

function [R_gl, ex, ey, ez] = ...
    build_local_frame(xyz1, xyz2, refVec)

    ex = xyz2 - xyz1;
    ex = ex / norm(ex);

    refVec = refVec(:).';
    refVec = refVec / norm(refVec);

    % 避免参考方向与梁轴线接近平行
    if abs(dot(ex, refVec)) > 0.95

        refVec = [0 1 0];

        if abs(dot(ex, refVec)) > 0.95
            refVec = [1 0 0];
        end
    end

    ez = cross(ex, refVec);
    ez = ez / norm(ez);

    ey = cross(ez, ex);
    ey = ey / norm(ey);

    ex = ex(:);
    ey = ey(:);
    ez = ez(:);

    % 局部矢量到全局矢量：
    % v_global = R_gl * v_local
    R_gl = [ex, ey, ez];

end


%% ========================================================================
% 局部函数：建立 12×12 坐标变换矩阵
%% ========================================================================

function T_gl = build_T_matrix(R_gl)

    T_node = blkdiag(R_gl, R_gl);

    T_gl = blkdiag(T_node, T_node);

end


%% ========================================================================
% 局部函数：一维 Gauss 积分
%% ========================================================================

function [xi, weight] = gauss_rule_1d(nGP)

    switch nGP

        case 1

            xi = 0;
            weight = 2;

        case 2

            xi = [-1; 1] / sqrt(3);
            weight = [1; 1];

        case 3

            xi = [-sqrt(3/5); 0; sqrt(3/5)];
            weight = [5/9; 8/9; 5/9];

        otherwise

            error('当前支持的 Gauss 积分点数为 1、2 或 3。');

    end

end