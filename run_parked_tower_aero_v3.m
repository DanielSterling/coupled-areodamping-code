function [model,response]=run_parked_tower_aero_v3(overrides)
% 例：run_parked_tower_aero_v3(struct('aeroChannels',6,'enableTimeResponse',true))
root=fileparts(mfilename('fullpath')); addpath(fullfile(root,'functions'),'-begin');
if nargin<1, overrides=struct(); end
cfg=parked_model_config(overrides); model=build_parked_model_v3(cfg); response=[];
fprintf('V3: %d aero channels; Ka=%d; mean=%d; wind=%.3f m/s\n', ...
 cfg.aeroChannels,cfg.enableAeroStiffness,cfg.includeMeanAeroLoad,cfg.inflowSpeed);
disp('Undamped structural frequencies (Hz):'); disp(model.modal.naturalFreqUndamped_Hz(1:cfg.numModesRetained).');
disp('Static tower-top displacement [m,m,m,rad,rad,rad]:'); disp(model.staticDisplacementTop.');
disp('  Gravity contribution:'); disp(model.staticDisplacementGravityTop.');
disp('  Mean-aero contribution (same effective stiffness):'); disp(model.staticDisplacementAeroTop.');
fprintf('max real state eigenvalue = %.9g /s\n',model.maxStateRealPart);
if model.maxStateRealPart>1e-8
 warning('model:Unstable','当前耦合线性系统存在增长模态；未裁剪负阻尼或强制稳定。');
end
if cfg.enableTimeResponse, response=solve_parked_response_v3(model); end
if cfg.saveResults
 if ~isfolder(cfg.resultFolder), mkdir(cfg.resultFolder); end
 save(fullfile(cfg.resultFolder,'parked_model_v3.mat'),'model','-v7.3');
 staticTable=table({'SS_X';'FA_Y';'Z';'thetaX';'thetaY';'thetaZ'}, ...
  {'m';'m';'m';'rad';'rad';'rad'},model.staticDisplacementAeroTop, ...
  model.staticDisplacementGravityTop,model.staticDisplacementTop, ...
  'VariableNames',{'DOF','Unit','MeanAero','Gravity','Total'});
 writetable(staticTable,fullfile(cfg.resultFolder,'tower_top_static_v3.csv'));
 writematrix(model.aero.C6,fullfile(cfg.resultFolder,'C6_aero.csv'));
 writematrix(model.aero.Jq6,fullfile(cfg.resultFolder,'Jq6_aero.csv'));
 writematrix(model.aero.Ka6,fullfile(cfg.resultFolder,'Ka6_aero.csv'));
 if ~isempty(response)
  save(fullfile(cfg.resultFolder,'parked_response_v3.mat'),'response','-v7.3');
  tableOut=array2table([response.time_s,response.towerTopPerturbation(:,[2,1]), ...
      response.towerTopTotal(:,[2,1])],'VariableNames', ...
      {'Time_s','FA_perturbation_m','SS_perturbation_m','FA_total_m','SS_total_m'});
  writetable(tableOut,fullfile(cfg.resultFolder,'tower_top_displacement_v3.csv'));
 end
end
if cfg.enablePlot && ~isempty(response)
 f=figure('Color','w'); tiledlayout(2,1);
 nexttile; plot(response.time_s,response.towerTopPerturbation(:,[2,1])); grid on;
 xlabel('Time (s)'); ylabel('Perturbation (m)'); legend('FA','SS');
 nexttile; plot(response.time_s,response.towerTopTotal(:,[2,1])); grid on;
 xlabel('Time (s)'); ylabel('Total incl. static offset (m)'); legend('FA','SS');
 if cfg.saveResults, exportgraphics(f,fullfile(cfg.resultFolder,'tower_top_response_v3.png'),'Resolution',180); end
end
end
