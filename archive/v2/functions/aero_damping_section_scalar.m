function out = aero_damping_section_scalar(Vx, Vy, c_chord, rho, Cl, Cd, dCl_dalpha, dCd_dalpha)
%AERO_DAMPING_SECTION_SCALAR  截面2×2气动阻尼矩阵（全标量输入版）
%
% 坐标说明：
%   本函数中的 Vx, Vy 不是全局坐标速度分量，
%   而是当前叶片叶根坐标系下截面平面 X_root-Y_root 内的速度分量。
%
%   Z_root 为叶片展向，当前二维翼型气动模型不使用 Vz_root。
%
%   out.Ca2 的自由度顺序为：
%       [X_root direction;
%        Y_root direction]
% 输入（全部为标量）：
%   Vx, Vy          : 风速在 X_B, Y_B 方向分量 (m/s)
%   c_chord         : 弦长 c (m)
%   rho             : 空气密度 rho (kg/m^3)
%   Cl, Cd          : 平衡点升力/阻力系数
%   dCl_dalpha      : ∂Cl/∂α (单位：1/rad)
%   dCd_dalpha      : ∂Cd/∂α (单位：1/rad)
%
% 输出 out：
%   out.cxx, out.cxy, out.cyx, out.cyy : 单位展向长度阻尼（N·s/m^2）
%   out.Ca2 (2x2)                      : 截面2×2阻尼矩阵
%   out.W                              : W = sqrt(Vx^2+Vy^2)
%
% 符号约定：
%   Ca2 = d(f_perLength)/d(V_relative)，空气速度固定时
%   delta_f_perLength = -Ca2 * v_sectionXY；装配到方程左侧为 +C*qdot。
%   f_perLength 单位 N/m；沿物理展向长度积分后得到总阻尼。
%   XY 坐标轴不随桨距/气动扭角旋转；这两个角只进入系数查询的攻角。
%
% 注意：
%   - 公式中有 1/W，若 W 很小会发散；做了下限保护。
%   - dCl_dalpha、dCd_dalpha 必须是“对弧度”的导数。

    % --- 相对风速大小 ---
    W = sqrt(Vx^2 + Vy^2);
    W_min = 1e-6;           % 防止除零
    W_eff = max(W, W_min);

    % --- 公共系数 0.5 * rho * c / W ---
    K = 0.5 * rho * c_chord / W_eff;

    Vx2  = Vx^2;
    Vy2  = Vy^2;
    VxVy = Vx * Vy;

    % ===== 按论文公式逐项实现 =====
    cxx = K * ( ...
          VxVy * dCd_dalpha ...
        + Vy2  * dCl_dalpha ...
        + (2*Vx2 + Vy2)*Cd ...
        + VxVy * Cl );

    cxy = K * ( ...
        - Vx2  * dCd_dalpha ...
        - VxVy * dCl_dalpha ...
        + (2*Vy2 + Vx2)*Cl ...
        + VxVy * Cd );

    cyx = K * ( ...
        - VxVy * dCl_dalpha ...
        + Vy2  * dCd_dalpha ...
        - (2*Vx2 + Vy2)*Cl ...
        + VxVy * Cd );

    cyy = K * ( ...
        - VxVy * dCd_dalpha ...
        + Vx2  * dCl_dalpha ...
        + (2*Vy2 + Vx2)*Cd ...
        - VxVy * Cl );

    % --- 输出 ---
    out = struct();
    out.cxx = cxx;
    out.cxy = cxy;
    out.cyx = cyx;
    out.cyy = cyy;

    out.W = W;

    out.Ca2 = [cxx, cxy;
               cyx, cyy];
end