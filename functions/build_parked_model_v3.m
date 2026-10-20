function out=build_parked_model_v3(cfg)
% 结构沿用原塔筒/RNA；新增参考气动力、位移切线及线性静平衡。
d=cfg.dataFolder;
[raw,~]=read_xlsx_numeric(fullfile(d,'tatong.xlsx'),1);
[nodes,elems,conn]=tower_elem_props_from_Tatong_gauss_FE(raw); top=nodes.xyz(end,:).';
rna=build_rna_lumped_mass_at_tower_top(fullfile(d,'nacelle.xlsx'), ...
 fullfile(d,'hub.xlsx'),fullfile(d,'blade_rna.xlsx'),top,cfg.bladeAzimuthSet_rad,cfg.bladePitchSet_rad);
if rna.mass<=0 || any(eig((rna.massMatrix6+rna.massMatrix6.')/2)<=0)
 error('model:RNA','RNA质量矩阵非正定，请核对质量、惯量和移轴参考点。');
end
model=assemble_tower_rna_structure(nodes,elems,conn,rna.massMatrix6);
free=(7:model.numDofTotal).'; K=model.K_global(free,free); M=model.M_global(free,free);
modal=solve_tower_rna_undamped_modes(K,M,cfg.numModesRequested); n=cfg.numModesRetained;
if n>modal.numModesSolved, error('model:ModeCount','保留模态数超过实际求解阶数。'); end
zeta=cfg.structuralDampingRatio(:);
if isscalar(zeta), zeta=repmat(zeta,n,1); end
if numel(zeta)~=n, error('model:DampingSize','结构阻尼比需为标量或保留模态数个元素。'); end
db=build_aerofoil_database(fullfile(d,'AerofoilDatabase.xlsx'));
sets=build_aerofoil_section_set(db,fullfile(d,'AerofoilSectionSet.xlsx'));
blade=read_blade_aero_geometry(fullfile(d,'blade_geo.xlsx'));
op.windVelocity_global=cfg.inflowSpeed*[cos(cfg.inflowAngle_rad);sin(cfg.inflowAngle_rad);0];
op.airDensity=cfg.airDensity; op.numGaussPoint=cfg.numGaussPoint;
aero=build_parked_aero_linearization(blade,rna,top,op,db,sets);
if cfg.aeroChannels==4, ids=[1,2,4,5]; else, ids=1:6; end
P=model.towerTopSelectionMatrix(:,free); Pa=P(ids,:);
Ct=aero.C6(ids,ids); Jt=aero.Jq6(ids,ids);
if ~cfg.enableAeroStiffness, Jt=zeros(size(Jt)); end
Q0=aero.Q0(ids); if ~cfg.includeMeanAeroLoad, Q0=zeros(size(Q0)); end
C=Pa.'*sparse(Ct)*Pa; Ka=-Pa.'*sparse(Jt)*Pa; Keff=K+Ka; F0=Pa.'*Q0;
% Gravity is a full structural load, independent of active aero channels.
gravity=build_parked_gravity_load(model,rna,cfg.gravityAcceleration_m_s2);
Fg=zeros(size(F0));
if cfg.includeGravity, Fg=gravity.forceGlobal(free); end
staticParts=Keff\[F0,Fg];
qstatic=sum(staticParts,2);
Fstatic=F0+Fg;
if any(~isfinite(qstatic)), error('model:Static','线性静平衡求解失败。'); end
staticResidual=norm(Keff*qstatic-Fstatic)/max(1,norm(Fstatic));
if staticResidual>1e-8, error('model:Static','线性静平衡残差过大。'); end
Phi=modal.modeShapeUndamped(:,1:n); H=P*Phi;
Mm=Phi.'*M*Phi; Km=Phi.'*Keff*Phi; Cam=Phi.'*C*Phi;
Csm=diag(2*zeta.*modal.naturalOmegaUndamped_rad_s(1:n).*diag(Mm));
A=[zeros(n),eye(n);-Mm\Km,-Mm\(Cam+Csm)]; B=[zeros(n,6);Mm\H.'];
out.config=cfg; out.structure=model; out.nodes=nodes; out.rna=rna;
out.freeDofId=free; out.K_free=K; out.M_free=M; out.modal=modal; out.Phi=Phi;
out.aero=aero; out.activeAeroDof=ids; out.topSelector=P;
out.C_free=C; out.Ka_free=Ka; out.Keffective_free=Keff;
out.Q0_active=Q0; out.Ctop_active=Ct; out.Jtop_active=Jt;
out.staticDisplacementFree=qstatic; out.staticDisplacementTop=P*qstatic;
out.staticResidual=staticResidual;
out.gravity=gravity; out.gravity.enabled=logical(cfg.includeGravity);
out.meanAeroForceFree=F0; out.gravityForceFree=Fg; out.staticForceFree=Fstatic;
out.staticDisplacementAeroFree=staticParts(:,1);
out.staticDisplacementGravityFree=staticParts(:,2);
out.staticDisplacementAeroTop=P*staticParts(:,1);
out.staticDisplacementGravityTop=P*staticParts(:,2);
out.M_modal=Mm; out.K_modal=Km; out.Ca_modal=Cam; out.Cs_modal=Csm;
out.A=A; out.B=B; out.H=H; out.stateEigenvalues=eig(A);
out.maxStateRealPart=max(real(out.stateEigenvalues));
% 仅作为未变形切线在所得静位移处的适用性诊断；不重新更新切线/质量。
probeQ=zeros(6,1); probeQ(ids)=out.staticDisplacementTop(ids);
probe=evaluate_parked_aero_force(probeQ,zeros(6,1),aero,db,sets);
predicted=aero.Q0(ids)+aero.Jq6(ids,ids)*probeQ(ids);
out.staticForceProbeRelativeError=norm(probe(ids)-predicted)/max(1,norm(probe(ids)));
out.referenceDescription=['Frozen undeformed tangent. Static offset solves the same affine model; ' ...
 'optional tower/RNA reference gravity load; no gravity load tangent, geometric stiffness, ' ...
 'nonlinear equilibrium iteration or deformed-state relinearization.'];
end
