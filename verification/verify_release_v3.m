function report=verify_release_v3
% 独立验证实际my260：系数斜率、速度/位移切线、静平衡和时程。
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'functions'),fullfile(root,'verification'),'-begin');
% Preserve the original gravity-free four-channel regression.
cfg=parked_model_config(struct('saveResults',false,'enablePlot',false, ...
 'aeroChannels',4,'includeGravity',false));
model=build_parked_model_v3(cfg);
db=build_aerofoil_database(fullfile(cfg.dataFolder,'AerofoilDatabase.xlsx'));
sets=build_aerofoil_section_set(db,fullfile(cfg.dataFolder,'AerofoilSectionSet.xlsx'));
aero=model.aero; z=zeros(6,1);
f=@(q,v)evaluate_parked_aero_force(q,v,aero,db,sets);
Q=f(z,z); report.meanForceRelativeError=rel(Q,aero.Q0);
assert(report.meanForceRelativeError<1e-12);
h=[1e-5;1e-5;1e-5;1e-7;1e-7;1e-7]; C=zeros(6); J=zeros(6); C2=C; J2=J;
for j=1:6
 e=z; e(j)=h(j);
 C(:,j)=-(f(z,e)-f(z,-e))/(2*h(j)); J(:,j)=(f(e,z)-f(-e,z))/(2*h(j));
 C2(:,j)=-(f(z,e/2)-f(z,-e/2))/h(j); J2(:,j)=(f(e/2,z)-f(-e/2,z))/h(j);
end
report.velocityTangentError=rel(aero.C6,C);
report.positionTangentError=rel(aero.Jq6,J);
report.velocityStepHalvingError=rel(C,C2); report.positionStepHalvingError=rel(J,J2);
assert(report.velocityTangentError<1e-7 && report.positionTangentError<1e-6);
assert(report.velocityStepHalvingError<1e-7 && report.positionStepHalvingError<1e-6);
assert(norm(aero.Jq6(:,1:3),'fro')==0);
ids=[1 2 4 5]; report.fourSixError=rel(aero.C6(ids,ids),aero.rotor4.dampingMatrix4);
assert(report.fourSixError<1e-12);
assert(rel(sum(aero.Q0_byBlade,2),aero.Q0)<1e-14);
assert(rel(sum(aero.Jq6_byBlade,3),aero.Jq6)<1e-14);
% 角点、区间内部、常值延拓及周期端点的明确约定。
aa=[0;1;2]; vv=[0;1;3]*[1,2,3];
[~,~,~,ds]=linear_polar_query(aa,vv,[-1,0,0.9,1,1.1,2,3],false);
assert(max(abs(ds-[0,0.5,1,1.5,2,1,0]))<1e-14);
[v,~,~,ds]=linear_polar_query([-pi;0;pi],[0;1;0]*[1,2,3],[-pi,pi,3*pi],true);
assert(max(abs(v))<1e-14 && max(abs(ds))<1e-14);
% 旧版的人工几何/光滑极曲线验证也必须继续通过。
verify_parked_aero_4dof;
% 非零整体倾斜、独立桨距及非水平风下验证新增位置切线。
rna=model.rna; rot=expm(skew([0.09;-0.04;0.02])); top=rna.towerTop_xyz;
for b=1:3
 g=rna.threeBlades.blade(b); pold=g.pitch_rad; pnew=deg2rad(70+7*b);
 g.rotRootToGlobal=rot*g.rotRootToGlobal*rz(pold)*rz(-pnew);
 g.rootReference_xyz=top+rot*(g.rootReference_xyz-top); g.pitch_rad=pnew;
 rna.threeBlades.blade(b)=g;
end
op=aero.operatingPoint; op.windVelocity_global=[-3;-5;0.7];
tilted=build_parked_aero_linearization(aero.blade,rna,top,op,db,sets);
ft=@(q,v)evaluate_parked_aero_force(q,v,tilted,db,sets);
jt=zeros(6); ct=zeros(6);
for j=1:6
 e=z; e(j)=h(j); jt(:,j)=(ft(e,z)-ft(-e,z))/(2*h(j)); ct(:,j)=-(ft(z,e)-ft(z,-e))/(2*h(j));
end
report.tiltedPositionError=rel(tilted.Jq6,jt); report.tiltedVelocityError=rel(tilted.C6,ct);
assert(report.tiltedPositionError<1e-6 && report.tiltedVelocityError<1e-7);
op.windVelocity_global=zeros(3,1);
zero=build_parked_aero_linearization(aero.blade,rna,top,op,db,sets);
assert(norm(zero.C6,'fro')==0 && norm(zero.Jq6,'fro')==0 && norm(zero.Q0)==0);
report.staticResidual=model.staticResidual;
assert(report.staticResidual<1e-8);
report.staticForceProbeRelativeError=model.staticForceProbeRelativeError;
report.staticDisplacementTop=model.staticDisplacementTop;
report.C4=aero.C6(ids,ids); report.Jq4=aero.Jq6(ids,ids); report.Q0=aero.Q0;
report.frequency_Hz=model.modal.naturalFreqUndamped_Hz(1:6);
% 位移反馈关闭后严格退化至结构刚度；六通道可正常建模。
cfgOff=cfg; cfgOff.enableAeroStiffness=false; cfgOff.includeMeanAeroLoad=false;
off=build_parked_model_v3(cfgOff);
assert(norm(off.Ka_free,'fro')==0 && norm(off.staticDisplacementFree)==0);
cfg6=cfg; cfg6.aeroChannels=6; six=build_parked_model_v3(cfg6);
assert(rel(six.Ctop_active,aero.C6)<1e-14 && six.staticResidual<1e-8);
report.maxRealEigenvalue_4ch=model.maxStateRealPart;
report.maxRealEigenvalue_6ch=six.maxStateRealPart;
% 载荷折点刻意不与0.02s输出网格对齐；与高精度ODE45独立对照。
tmp=[tempname,'.txt']; fid=fopen(tmp,'w');
fprintf(fid,'NUMPOINTS 4\n0 0\n0.033 1000\n0.073 -300\n0.111 0\n'); fclose(fid);
cleanup=onCleanup(@()delete(tmp));
model.config.excitationFile=tmp; model.config.totalTime_s=0.4;
r=solve_parked_response_v3(model);
impact=read_impact_load(tmp,0);
ode=@(t,x)model.A*x+model.B*cfg.excitationDirection6*interp1(impact.timeGlobal,impact.force,t,'linear',0);
[~,xo]=ode45(ode,r.time_s,zeros(12,1),odeset('RelTol',1e-10,'AbsTol',1e-13,'MaxStep',0.001));
report.timeResponseRelativeError=norm(r.stateHistory-xo,'fro')/max(1e-14,norm(xo,'fro'));
assert(report.timeResponseRelativeError<1e-6);
assert(rel(r.towerTopTotal-r.towerTopPerturbation,repmat(model.staticDisplacementTop.',numel(r.time_s),1))<1e-14);
% 原外部载荷若可读，做完整400秒实机回归。
report.actualImpactRun=false;
if isfile(cfg.excitationFile)
 cfgRun=cfg; cfgRun.enableTimeResponse=true;
 try
  read_impact_load(cfgRun.excitationFile,0);
 catch err
  error('release:Impact','随包Impact不可用：%s',err.message);
 end
 model.config=cfgRun; realResponse=solve_parked_response_v3(model);
 report.actualImpactRun=true; report.actualTimePoints=numel(realResponse.time_s);
 report.actualPeakPerturbationFA_SS=max(abs(realResponse.towerTopPerturbation(:,[2,1])),[],1);
end
report.status='PASS'; report.matlabVersion=version; report.date='2026-10-16';
fid=fopen(fullfile(root,'verification','verification_results.json'),'w','n','UTF-8');
fprintf(fid,'%s\n',jsonencode(report,'PrettyPrint',true)); fclose(fid);
disp(jsonencode(report));
end
function v=rel(a,b)
v=norm(a(:)-b(:))/max(1,norm(b(:)));
end
function S=skew(a)
S=[0,-a(3),a(2);a(3),0,-a(1);-a(2),a(1),0];
end
function R=rz(a)
R=[cos(a),-sin(a),0;sin(a),cos(a),0;0,0,1];
end
