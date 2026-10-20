function response=solve_parked_response_v3(model)
cfg=model.config; impact=read_impact_load(cfg.excitationFile,cfg.excitationStartTime_s);
t=(0:cfg.timeStep_s:cfg.totalTime_s).';
if t(end)<cfg.totalTime_s, t(end+1)=cfg.totalTime_s; end
knots=impact.timeGlobal(impact.timeGlobal>=0 & impact.timeGlobal<=cfg.totalTime_s);
ti=unique([t;knots]); [~,oi]=ismember(t,ti);
F=interp1(impact.timeGlobal,impact.force,ti,'linear',0)*cfg.excitationDirection6(:).';
n=cfg.numModesRetained; ns=2*n;
aug=[model.A,model.B,zeros(ns,6);zeros(6,ns+6),eye(6);zeros(6,ns+12)];
x=zeros(numel(ti),ns); cache=containers.Map('KeyType','char','ValueType','any');
for k=1:numel(ti)-1
 dt=ti(k+1)-ti(k); key=sprintf('%.17g',dt);
 if isKey(cache,key), E=cache(key); else, E=expm(aug*dt); cache(key)=E; end
 slope=(F(k+1,:)-F(k,:)).'/dt;
 x(k+1,:)=(E(1:ns,1:ns)*x(k,:).'+E(1:ns,ns+(1:6))*F(k,:).' ...
    +E(1:ns,ns+6+(1:6))*slope).';
end
x=x(oi,:); F=F(oi,:); eta=x(:,1:n); etad=x(:,n+1:end);
if any(~isfinite(x(:))), error('model:ResponseOverflow','响应出现非有限值，请检查系统不稳定性。'); end
top=eta*model.H.'; velocity=etad*model.H.';
acc=(model.M_modal\((F*model.H).'-(model.Ca_modal+model.Cs_modal)*etad.'-model.K_modal*eta.')).';
response.time_s=t; response.impact=impact; response.forceTop=F;
response.stateHistory=x; response.modalAcceleration=acc;
response.towerTopPerturbation=top;
response.towerTopTotal=top+model.staticDisplacementTop.';
response.towerTopVelocity=velocity; response.towerTopAcceleration=acc*model.H.';
ids=model.activeAeroDof;
response.aeroLoadPerturbation=top(:,ids)*model.Jtop_active.'-velocity(:,ids)*model.Ctop_active.';
response.aeroLoadMean=model.Q0_active+model.Jtop_active*model.staticDisplacementTop(ids);
response.aeroLoadTotal=response.aeroLoadPerturbation+response.aeroLoadMean.';
response.displacementDescription='Perturbation about affine-model static offset; total includes offset.';
end
