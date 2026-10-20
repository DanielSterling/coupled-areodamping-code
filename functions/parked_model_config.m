function cfg=parked_model_config(overrides)
% 可通过函数参数覆盖设置，无需编辑/保存源文件。
root=fileparts(fileparts(mfilename('fullpath')));
cfg.dataFolder=fullfile(root,'my260');
cfg.aeroChannels=6; % 4保留原通道近似；6启用全部塔顶气动通道。
cfg.enableAeroStiffness=true;
cfg.includeMeanAeroLoad=true;
cfg.includeGravity=true; % Tower self-weight and eccentric RNA weight.
cfg.gravityAcceleration_m_s2=9.81;
cfg.numModesRequested=16; cfg.numModesRetained=6;
cfg.structuralDampingRatio=0;
cfg.bladeAzimuthSet_rad=deg2rad([0;120;240]);
cfg.bladePitchSet_rad=deg2rad([90;90;90]);
cfg.inflowSpeed=6; cfg.inflowAngle_rad=deg2rad(240);
cfg.airDensity=1.29; cfg.numGaussPoint=3;
cfg.enableTimeResponse=false; % 开箱先算模型；显式开启时程，要求可读Impact文件。
cfg.excitationFile=fullfile(root,'inputs','FA_modal_excitation.txt');
cfg.excitationStartTime_s=0; cfg.excitationDirection6=[0;-1;0;0;0;0];
cfg.timeStep_s=0.05; cfg.totalTime_s=1000;
cfg.enablePlot=true; cfg.saveResults=true; cfg.resultFolder='';
if nargin>0
    fields=fieldnames(overrides);
    for k=1:numel(fields)
        if ~isfield(cfg,fields{k}), error('model:Config','未知设置%s。',fields{k}); end
        cfg.(fields{k})=overrides.(fields{k});
    end
end
if ~ismember(cfg.aeroChannels,[4,6]), error('model:Config','aeroChannels只能为4或6。'); end
validateattributes(cfg.numModesRequested,{'numeric'},{'scalar','integer','positive','finite'});
validateattributes(cfg.numModesRetained,{'numeric'},{'scalar','integer','positive','finite'});
validateattributes(cfg.inflowSpeed,{'numeric'},{'scalar','nonnegative','finite'});
validateattributes(cfg.inflowAngle_rad,{'numeric'},{'scalar','real','finite'});
validateattributes(cfg.airDensity,{'numeric'},{'scalar','positive','finite'});
validateattributes(cfg.gravityAcceleration_m_s2,{'numeric'},{'scalar','real','nonnegative','finite'});
validateattributes(cfg.includeGravity,{'logical','numeric'},{'scalar','real','finite','binary'});
validateattributes(cfg.timeStep_s,{'numeric'},{'scalar','positive','finite'});
validateattributes(cfg.totalTime_s,{'numeric'},{'scalar','positive','finite'});
validateattributes(cfg.excitationStartTime_s,{'numeric'},{'scalar','real','finite'});
validateattributes(cfg.bladePitchSet_rad,{'numeric'},{'vector','numel',3,'real','finite'});
validateattributes(cfg.bladeAzimuthSet_rad,{'numeric'},{'vector','numel',3,'real','finite'});
validateattributes(cfg.excitationDirection6,{'numeric'},{'vector','numel',6,'real','finite'});
validateattributes(cfg.structuralDampingRatio,{'numeric'},{'vector','real','finite','nonnegative'});
if ~ismember(cfg.numGaussPoint,[2,3]), error('model:Config','高斯点数为2或3。'); end
if isempty(cfg.resultFolder)
    cfg.resultFolder=fullfile(root,'results',sprintf('v3_%dch_Ka%d_mean%d_g%d', ...
        cfg.aeroChannels,cfg.enableAeroStiffness,cfg.includeMeanAeroLoad,cfg.includeGravity));
end
end
