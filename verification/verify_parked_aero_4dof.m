function verify_parked_aero_4dof
%VERIFY_PARKED_AERO_4DOF 不依赖机型Excel的物理验证，使用人工光滑极曲线。
% 在MATLAB中运行本文件；验证实际build_parked_rotor_aero_damping_4dof。
% 人工几何/极曲线仅用于检查算法，不能解释为my260的阻尼结果。
projectFolder = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectFolder, 'functions'));

%% 1. 构造停机三叶片：倾角、锥角、偏置和独立桨距均非零
bladeData.distanceAlongBlade = [0; 30; 70; 120];
bladeData.spanPositionFromRoot = [0; 30; 70; 120];
bladeData.chord = [5; 4; 2; 1];
bladeData.aeroTwist_rad = deg2rad([10; 5; 1; 0]);
bladeData.thickness_pct = 30*ones(4, 1);
bladeData.foilSection = ones(4, 1);
towerTop = [0; 0; 150];
hubReference = towerTop + [1; 6; 3];
tilt = deg2rad(6); cone = deg2rad(3);
shaft = [0; cos(tilt); sin(tilt)];
ex = -shaft;
ez = [0; 0; 1] - ex*ex(3); ez = ez/norm(ez);
R0 = [ex, cross(ez, ex), ez];
pitchSet = deg2rad([90; 78; 85]);
azimuthSet = deg2rad([0; 120; 240]);
for bladeIdx = 1:3
    referenceR = R0*rotation_x(azimuthSet(bladeIdx))*rotation_y(cone);
    rna.threeBlades.blade(bladeIdx).pitch_rad = pitchSet(bladeIdx);
    rna.threeBlades.blade(bladeIdx).azimuth_rad = azimuthSet(bladeIdx);
    rna.threeBlades.blade(bladeIdx).rotRootToGlobal = ...
        referenceR*rotation_z(-pitchSet(bladeIdx));
    rna.threeBlades.blade(bladeIdx).rootReference_xyz = ...
        hubReference + referenceR*[0; 0; 3];
end
op.windVelocity_global = 20*[cosd(240); sind(240); 0];
op.airDensity = 1.29;
op.numGaussPoint = 3;
database.query = @(name, alpha) mock_polar(alpha);
sectionSet.sectionNumber = 1;
sectionSet.foilNameList = "mock";
sectionSet.tcPercentList = 30;
out = build_parked_rotor_aero_damping_4dof( ...
    bladeData, rna, towerTop, op, database, sectionSet);
C4 = out.dampingMatrix4;
assert_close(C4, sum(out.dampingMatrix4_byBlade, 3), 1e-13, '逐叶片求和');

%% 2. 用六自由度刚体映射独立组装：直接C4应等于C6([1 2 4 5],[1 2 4 5])
% 这里只在验证中计算C6，生产函数没有生成C6。
C6 = zeros(6);
q4dot = [0.2; -0.1; 0.003; -0.002];
sectionPower = 0;
for bladeIdx = 1:3
    blade = out.bladeAeroOut(bladeIdx);
    T = blade.rotationReferenceToGlobal(:, 1:2);
    gp = blade.gaussPointInfo;
    for k = 1:numel(gp.integrationWeight_m)
        r = gp.rTowerTopToSection_global(:, k);
        B6 = T.'*[eye(3), -skew_matrix(r)];
        c2 = gp.sectionDampingMatrix2(:, :, k);
        w = gp.integrationWeight_m(k);
        C6 = C6 + B6.'*c2*B6*w;
        v = gp.sectionVelocityMap4(:, :, k)*q4dot;
        sectionPower = sectionPower + v.'*c2*v*w;
        velocityExpected = T.'*([q4dot(1:2); 0] ...
            + cross([q4dot(3:4); 0], r));
        assert_close(v, velocityExpected, 1e-12, '截面速度叉乘');
    end
end
index4 = [1; 2; 4; 5];
assert_close(C4, C6(index4, index4), 1e-12, '4/6自由度矩阵关系');
assert_close(q4dot.'*C4*q4dot, sectionPower, 1e-12, '功率等效');

%% 3. 对四个塔顶速度施加正负扰动，从非线性气动力求导
% 平移步长m/s，转动步长rad/s；逐列中心差分。
step = [1e-5; 1e-5; 1e-7; 1e-7];
forceDerivative = zeros(4);
for j = 1:4
    dq = zeros(4, 1); dq(j) = step(j);
    forceDerivative(:, j) = (tower_force4(dq, out, bladeData, op) ...
        - tower_force4(-dq, out, bladeData, op))/(2*step(j));
end
assert_close(-forceDerivative, C4, 1e-7, '非线性载荷导数=-C4');

%% 4. 两点/三点Gauss对当前分段线性插值模型应一致
op.numGaussPoint = 2;
out2 = build_parked_rotor_aero_damping_4dof( ...
    bladeData, rna, towerTop, op, database, sectionSet);
assert_close(out2.dampingMatrix4, C4, 1e-12, '2/3点Gauss');

%% 5. 零风速下截面气动力速度雅可比为零
op.windVelocity_global = zeros(3, 1);
outZero = build_parked_rotor_aero_damping_4dof( ...
    bladeData, rna, towerTop, op, database, sectionSet);
assert_close(outZero.dampingMatrix4, zeros(4), 1e-13, '零风速');

disp('人工极曲线/几何得到的示例C4，非my260机型结果：');
disp(C4);
disp('全部四自由度气动阻尼物理验证通过。');
end

function Q4 = tower_force4(q4dot, out, bladeData, op)
Q6 = zeros(6, 1);
for b = 1:numel(out.bladeAeroOut)
    blade = out.bladeAeroOut(b);
    T = blade.rotationReferenceToGlobal(:, 1:2);
    V0 = blade.windVelocity_reference(1:2);
    gp = blade.gaussPointInfo;
    for k = 1:numel(gp.integrationWeight_m)
        i = gp.segmentIndex(k);
        fraction = (gp.distanceAlongBlade(k)-bladeData.distanceAlongBlade(i)) ...
            /(bladeData.distanceAlongBlade(i+1)-bladeData.distanceAlongBlade(i));
        r = gp.rTowerTopToSection_global(:, k);
        % 独立使用cross求截面速度，不复用生产代码中的B4。
        vGlobal = [q4dot(1:2); 0] + cross([q4dot(3:4); 0], r);
        V = V0 - T.'*vGlobal;
        fi = section_force(V, blade.bladePitch_rad, ...
            bladeData.aeroTwist_rad(i), bladeData.chord(i), op.airDensity);
        fj = section_force(V, blade.bladePitch_rad, ...
            bladeData.aeroTwist_rad(i+1), bladeData.chord(i+1), op.airDensity);
        fGlobal = T*((1-fraction)*fi + fraction*fj);
        Q6 = Q6 + [fGlobal; cross(r, fGlobal)]*gp.integrationWeight_m(k);
    end
end
Q4 = Q6([1; 2; 4; 5]);
end

function f = section_force(V, pitch, twist, chord, rho)
alpha = atan2(V(1), V(2))-pitch-twist;
[Cl, Cd] = mock_polar(alpha);
f = 0.5*rho*chord*norm(V) ...
    * [Cd*V(1)+Cl*V(2); Cd*V(2)-Cl*V(1)];
end

function [Cl, Cd, Cm, dCl, dCd, dCm] = mock_polar(alpha)
Cl = 0.7*sin(alpha); Cd = 0.3+0.1*cos(2*alpha); Cm = zeros(size(alpha));
dCl = 0.7*cos(alpha); dCd = -0.2*sin(2*alpha); dCm = zeros(size(alpha));
end

function assert_close(actual, expected, tolerance, name)
relativeError = norm(actual(:)-expected(:))/max(1, norm(expected(:)));
assert(relativeError < tolerance, '%s失败：相对误差%.3e。', name, relativeError);
end

function S = skew_matrix(r)
S = [0 -r(3) r(2); r(3) 0 -r(1); -r(2) r(1) 0];
end

function R = rotation_x(a)
R = [1 0 0; 0 cos(a) -sin(a); 0 sin(a) cos(a)];
end

function R = rotation_y(a)
R = [cos(a) 0 sin(a); 0 1 0; -sin(a) 0 cos(a)];
end

function R = rotation_z(a)
R = [cos(a) -sin(a) 0; sin(a) cos(a) 0; 0 0 1];
end
