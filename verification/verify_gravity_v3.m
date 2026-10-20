function report=verify_gravity_v3
% Independent balance, static superposition, channel and response checks.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'functions'),'-begin');
cfg=parked_model_config(struct('saveResults',false,'enablePlot',false,'includeGravity',true));
m=build_parked_model_v3(cfg); g=cfg.gravityAcceleration_m_s2;
% Integrate the linearly interpolated element masses independently of M.
expectedMass=0;
% Recreate station properties rather than infer weight from the same matrix.
[raw,~]=read_xlsx_numeric(fullfile(cfg.dataFolder,'tatong.xlsx'),1);
[nodes,elems,~]=tower_elem_props_from_Tatong_gauss_FE(raw);
expectedMoment=zeros(3,1);
for e=1:numel(elems)
 el=elems(e); r1=el.xyz1(:); r2=el.xyz2(:); L=norm(r2-r1);
 expectedMass=expectedMass+L*(el.m1+el.m2)/2;
 for xi=[-1,1]/sqrt(3)
  N=[(1-xi)/2,(1+xi)/2]; r=r1*N(1)+r2*N(2);
  dm=(N*[el.m1;el.m2])*L/2;
  expectedMoment=expectedMoment+cross(r,[0;0;-g*dm]);
 end
end
ft=reshape(m.gravity.forceTowerGlobal,6,[]);
actualMoment=sum(ft(4:6,:)+cross(nodes.xyz.',ft(1:3,:),1),2);
report.towerForceError=relative(sum(ft(1:3,:),2),[0;0;-expectedMass*g]);
report.towerMomentError=relative(actualMoment,expectedMoment);
assert(report.towerForceError<1e-10 && report.towerMomentError<1e-10);
report.rnaReferenceRelativeError=m.gravity.rnaReferenceRelativeError;
assert(report.rnaReferenceRelativeError<1e-10);
report.staticResidual=m.staticResidual; assert(report.staticResidual<1e-8);
assert(relative(m.staticDisplacementTop, ...
 m.staticDisplacementAeroTop+m.staticDisplacementGravityTop)<1e-12);
cfg.includeGravity=false; off=build_parked_model_v3(cfg);
assert(norm(off.gravityForceFree)==0);
assert(relative(off.staticDisplacementTop,m.staticDisplacementAeroTop)<1e-12);
assert(relative(off.A,m.A)<1e-14 && relative(off.B,m.B)<1e-14);
cfg.includeGravity=true; cfg.aeroChannels=4; four=build_parked_model_v3(cfg);
assert(relative(four.gravityForceFree,m.gravityForceFree)<1e-14);
% Gravity must remain present when mean aerodynamic loading is disabled.
cfg.includeMeanAeroLoad=false; onlyGravity=build_parked_model_v3(cfg);
assert(norm(onlyGravity.meanAeroForceFree)==0);
assert(norm(onlyGravity.staticDisplacementGravityTop)>0);
assert(relative(onlyGravity.staticDisplacementTop,four.staticDisplacementGravityTop)<1e-12);
% Synthetic short impact: same perturbation, different constant offset.
tmp=[tempname,'.txt']; fid=fopen(tmp,'w');
assert(fid>=0); fprintf(fid,'NUMPOINTS 4\n0 0\n0.033 1000\n0.073 -300\n0.111 0\n'); fclose(fid);
cleanup=onCleanup(@()delete(tmp)); %#ok<NASGU>
m.config.excitationFile=tmp; m.config.totalTime_s=0.4;
off.config.excitationFile=tmp; off.config.totalTime_s=0.4;
r=solve_parked_response_v3(m); ro=solve_parked_response_v3(off);
assert(norm(r.towerTopPerturbation(1,:))==0);
assert(relative(r.towerTopTotal(1,:).',m.staticDisplacementTop)<1e-14);
assert(relative(r.towerTopPerturbation,ro.towerTopPerturbation)<1e-14);
assert(relative(r.towerTopTotal-ro.towerTopTotal, ...
 repmat(m.staticDisplacementGravityTop.',numel(r.time_s),1))<1e-12);
% Zero gravity acceleration is another exact regression path.
zero=build_parked_gravity_load(m.structure,m.rna,0);
assert(norm(zero.forceGlobal)==0);
report.staticDisplacementTop=m.staticDisplacementTop;
report.staticGravityTop=m.staticDisplacementGravityTop;
report.staticAeroTop=m.staticDisplacementAeroTop;
report.staticForceProbeRelativeError=m.staticForceProbeRelativeError;
report.status='PASS'; report.matlabVersion=version;
fid=fopen(fullfile(root,'verification','gravity_verification_results.json'),'w','n','UTF-8');
assert(fid>=0); fprintf(fid,'%s\n',jsonencode(report,'PrettyPrint',true)); fclose(fid);
disp(jsonencode(report));
end
function e=relative(a,b)
e=norm(a(:)-b(:))/max(1,norm(b(:)));
end
