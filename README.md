# coupled-areodamping-code

用于描述耦合气动反馈并预测塔筒顶部位移时程响应。当前主版本：my260 停机气动反馈模型 V3。

发布整理日期：2026-10-20；源码快照包含截至2026-10-19的本地修改。适用 MATLAB R2024b（本次验证环境）。

本包独立于旧工程，保留真实 my260 输入与原代码备份。建议完整解压到新目录，切换 MATLAB 当前目录到本包根目录后运行。不要把 reference_original 加入 MATLAB 路径，也不要对整个包使用 addpath(genpath(...))。

## 运行方式：无需编辑源文件

默认计算结构、气动切线、线性静平衡与模态状态矩阵，开启六通道位移反馈（含重力静平衡）：

```matlab
[model,response] = run_parked_tower_aero_v3;
```

同样可以运行 `main_parked_tower_aero_4dof_v3.m`。默认不计算时程，response 为空。

使用随包原始 Impact 载荷计算默认1000秒时程：

```matlab
[model,response] = run_parked_tower_aero_v3(struct('enableTimeResponse',true));
```

与原“纯阻尼、零平均载荷”假设对照（翼型导数仍采用修正版）：

```matlab
[model,response] = run_parked_tower_aero_v3(struct( ...
    'enableAeroStiffness',false,'includeMeanAeroLoad',false,'includeGravity',false, ...
    'enableTimeResponse',true));
```

启用六个塔顶气动通道：

```matlab
[model,response] = run_parked_tower_aero_v3(struct( ...
    'aeroChannels',6,'enableTimeResponse',true));
```

指定自己的载荷路径，也通过参数传入，不需要保存代码：

```matlab
run_parked_tower_aero_v3(struct('enableTimeResponse',true, ...
    'excitationFile','D:\your_folder\your_impact.txt'));
```

全部可覆盖设置见 `functions/parked_model_config.m`。改变风速、模态数等参数后如需保留多组结果，请同时指定不同 `resultFolder`；相同输出目录的同名文件会更新。

## 默认机型及工况

- 数据目录：my260，保留9个原始Excel；V3使用其中7个，未使用 blade_mass_stiff.xlsx 和 blade_pointmass.xlsx。
- 全局X=SS，Y=FA，Z向上；每节点自由度为 [uX,uY,uZ,thetaX,thetaY,thetaZ]。
- RNA刚性、塔筒柔性、固定塔底；RNA质量为737365 kg。
- 风速6 m/s，空气速度方向角240°，空气密度1.29 kg/m³；不是气象学“风从哪里来”的方位角。
- 方位角[0,120,240]°，桨距[90,90,90]°；倾角/锥角从hub.xlsx读取，当前均为0°。
- 结构求解16阶，保留前6阶动态模态，结构阻尼比0；不对负气动阻尼或增长模态做截断。
- 默认六通道；aeroChannels=4可选择[1,2,4,5]。四通道不约束结构uZ/thetaZ，只省略相应气动载荷与反馈。
- 时程输出间隔0.05 s，总长1000 s；对比Bladed时显式设totalTime_s=400；所有Impact折点另行加入内部积分网格。

## 新模型的精确定义

在未变形、零运动速度的参考状态建立一次切线：

`Qa(q,v) ≈ Q0 + Jq*q - Ca*v`，`Ka = -Jq`。

气动切线使用全局物理力及关于随动塔顶的全局物理力矩，转角为小转角增量。塔筒结构仍是小变形线性结构。

均匀风中纯平移不改变参考气动力，Jq前三列为零；转角改变截面坐标轴、风速投影、力方向及力臂，Jq后三列一般不为零。

选定气动通道后，通过塔顶选择矩阵P装配：

`Cfree = P'*Ca*P`，`Kafree = -P'*Jq*P`，`F0free = P'*Q0`。

线性静平衡：`(Ks + Kafree)*qstatic = F0free + Fgravity_free`（includeGravity=false 时重力项为零）。

动态扰动：`q = qstatic + Phi*eta`，满足

`Mm*eta_ddot + (Cs_m+Ca_m)*eta_dot + (Ks_m+Ka_m)*eta = H'*Fimpact`。

初始**扰动**位移/速度为零；初始总位移包含平均气动力和已启用重力产生的qstatic。

这是在同一未变形切线下求得的线性/仿射静平衡。它没有迭代非线性静平衡，没有按静变形重新线性化，也没有更新RNA质量矩阵。`staticForceProbeRelativeError`用重新计算姿态的气动力检查该近似，不能代替非线性静平衡求解。静平衡默认包含叶片平均气动力、塔筒自重和 RNA 偏心重力；不包含塔筒气动力、重力随姿态变化的切线或几何刚度。

如果 includeMeanAeroLoad=false，则参考平均气动力被移除；重力仍由includeGravity独立控制，气动刚度由enableAeroStiffness独立控制。

## 输出含义

各设置写入独立的 `results/v3_6ch_Ka1_mean1_g1` 等目录。

| 输出 | 内容 |
| --- | --- |
| parked_model_v3.mat | model结构体：结构矩阵、模态、全六通道气动切线、实际启用通道、静位移、状态矩阵及特征值 |
| C6_aero.csv | 全六通道参考气动阻尼；选择四通道运行时仍导出全矩阵供诊断 |
| Jq6_aero.csv | 全六通道气动力对位移的导数，行力/力矩，列位移/转角 |
| Ka6_aero.csv | `-Jq6`，移到动力学方程左侧的气动刚度 |
| parked_response_v3.mat | response结构体：扰动/总位移、速度、加速度、外载和气动载荷 |
| tower_top_displacement_v3.csv | 时间、FA扰动位移、SS扰动位移、FA总位移、SS总位移 |
| tower_top_response_v3.png | 扰动位移和总位移分开绘制 |

六通道矩阵行顺序[FX,FY,FZ,MX,MY,MZ]，列顺序[uX,uY,uZ,thetaX,thetaY,thetaZ]（Ca的列为相应速度）。

Ca分块单位依次为N·s/m、N·s、N·s、N·m·s。Jq/Ka分块单位为N/m、N、N、N·m，rad视为无量纲。任何矩阵都不代表阻尼比。

`response.aeroLoadPerturbation`为选定通道的`Jq*delta_q-Ca*delta_v`；`aeroLoadTotal`再加入线性静平衡点的气动载荷。该载荷不应再次添加到已经包含Ka/Ca的状态方程中，否则会重复计入。

## 验证与适用范围

```matlab
addpath('verification');
report = verify_release_v3;
```

验证结果写入 `verification/verification_results.json`。具体已执行结果见更新文档。

仍保留准定常二维极曲线、直线刚性叶片及中心线作用点假设。未新增叶片柔性、Cm反馈、非定常失速、风切变、湍流、塔筒涡激力或结构几何非线性。代码验证证明实现与声明的模型一致，不证明与Bladed或实测完全一致。

## 发布内容与版本对照

当前根目录为V3，`archive/v2/`保留本次整理时聊天工作目录中的V2主程序、通用函数和历史说明。V2说明中的默认值和验证环境为历史记录，以实际源码为准。

V2的文件载荷入口仍含原机器绝对路径，运行前须自行指定；V3默认使用随仓库提供的`inputs/FA_modal_excitation.txt`。不要将V2目录与V3同时加入MATLAB搜索路径。V2若复用根目录机型输入，请将`my260`放到V2入口旁，或显式修改其`turbineDataFolder`。

仓库包含9个机型Excel输入、Impact载荷、验证脚本及历史验证JSON；不包含计算结果、临时文件或修改前备份。`SHA256SUMS.txt`用于校验本次发布文件。


## 2026-10-17 重力静平衡更新

当前实际默认配置为6个气动通道，includeGravity=true，gravityAcceleration_m_s2=9.81。
重力用一致质量矩阵乘全体节点的竖向刚体加速度场生成：Fg=M*a_g。
塔筒和RNA贡献分别保存；RNA结果与[m*g; r_G cross (m*g)]核对，避免重复移轴。
重力保留全部结构自由度，不受4/6气动通道选择限制；静平衡使用全部自由自由度求解。
staticDisplacementAeroTop与staticDisplacementGravityTop为同一有效刚度下的线性分解，二者之和为总静位移。
gravity.forceGlobal等字段保存候选重力载荷，gravity.enabled和gravityForceFree反映是否实际施加。
结果目录新增_g0/_g1后缀；新增tower_top_static_v3.csv，包含六个自由度的单位和静位移分量。

```matlab
% 共同400秒时段：开启重力并计算响应
[model,response]=run_parked_tower_aero_v3(struct( ...
 'includeGravity',true,'enableTimeResponse',true,'totalTime_s',400));
% 回归原先只有平均气动力的静平衡
[oldModel,~]=run_parked_tower_aero_v3(struct('includeGravity',false));
% 验证重力合力、合矩、静平衡分解与动态初始状态
addpath('verification'); verify_gravity_v3;
```

这是未变形参考状态下的恒定重力载荷更新。未加入P-Delta、RNA重力姿态切线或变形后重新线性化；
因此本次更新改变静位移和总位移参考，动态扰动、频率及增长率保持原模型计算方式。
风向240度与Bladed风向30度的映射、位移轴符号及Impact载荷一致性仍需核对。


## 本次发布验证

本次MATLAB验证记录见`RELEASE_NOTES.md`。历史JSON文件仅代表对应验证脚本的工况，不等同于最新默认配置全部已验证。
