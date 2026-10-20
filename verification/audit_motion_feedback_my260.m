function report = audit_motion_feedback_my260(mode)
% Read-only model audit. Rebuilds the current 6 m/s, 240 deg parked case.
% Independent force evaluation uses global cross products and rotated geometry.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'functions'));
d = fullfile(root,'my260');
[raw,~] = read_xlsx_numeric(fullfile(d,'tatong.xlsx'),1);
[nodes,elems,conn] = tower_elem_props_from_Tatong_gauss_FE(raw);
top = nodes.xyz(end,:).';
rna = build_rna_lumped_mass_at_tower_top(fullfile(d,'nacelle.xlsx'), ...
    fullfile(d,'hub.xlsx'),fullfile(d,'blade_rna.xlsx'),top, ...
    deg2rad([0;120;240]),deg2rad([90;90;90]));
blade = read_blade_aero_geometry(fullfile(d,'blade_geo.xlsx'));
db = build_aerofoil_database(fullfile(d,'AerofoilDatabase.xlsx'));
sections = build_aerofoil_section_set(db,fullfile(d,'AerofoilSectionSet.xlsx'));
op.windVelocity_global = 6*[cosd(240);sind(240);0];
op.airDensity = 1.29; op.numGaussPoint = 3;
rotor = build_parked_rotor_aero_damping_4dof(blade,rna,top,op,db,sections);
% Isolate the polar derivative approximation using a much smaller angular step.
Csmall=zeros(4); derivativeMaxDifference=0; affected=0;
for b=1:3
 g=rotor.bladeAeroOut(b); alpha=g.sectionAoA_rad; da=1e-7;
 [lp,dp]=query_aerofoil_section(db,sections,blade.foilSection,blade.thickness_pct,alpha+da);
 [lm,dm]=query_aerofoil_section(db,sections,blade.foilSection,blade.thickness_pct,alpha-da);
 dl=(lp-lm)/(2*da); dd=(dp-dm)/(2*da);
 err=max(abs([dl-g.sectionLiftCoeff_dAlpha,dd-g.sectionDragCoeff_dAlpha]),[],2);
 derivativeMaxDifference=max(derivativeMaxDifference,max(err));
 affected=affected+sum(err>1e-6);
 cs=zeros(2,2,numel(alpha)); V=g.windVelocity_reference;
 for i=1:numel(alpha)
  s=aero_damping_section_scalar(V(1),V(2),blade.chord(i),op.airDensity, ...
      g.sectionLiftCoeff(i),g.sectionDragCoeff(i),dl(i),dd(i));
  cs(:,:,i)=s.Ca2;
 end
 gp=g.gaussPointInfo;
 for k=1:numel(gp.integrationWeight_m)
  i=gp.segmentIndex(k);
  t=(gp.distanceAlongBlade(k)-blade.distanceAlongBlade(i))/(blade.distanceAlongBlade(i+1)-blade.distanceAlongBlade(i));
  B=gp.sectionVelocityMap4(:,:,k);
  Csmall=Csmall+B.'*((1-t)*cs(:,:,i)+t*cs(:,:,i+1))*B*gp.integrationWeight_m(k);
 end
end
report.C4_smallAngleDerivative=Csmall;
report.polarDerivativeMaxDifference_perRad=derivativeMaxDifference;
report.polarDerivativeAffectedSections=affected;
if nargin>0 && strcmp(mode,'derivatives')
 disp(jsonencode(report)); return
end
z = zeros(6,1); idx = [1 2 4 5];
force = @(q,v) independent_force(q,v,rotor,blade,op,db,sections);
Q0 = force(z,z);
h = [1e-5;1e-5;1e-5;1e-7;1e-7;1e-7];
C6 = zeros(6); C6half = zeros(6); Jq6 = zeros(6); Jq6half = zeros(6);
for j = 1:6
    e=z; e(j)=h(j);
    C6(:,j)=-(force(z,e)-force(z,-e))/(2*h(j));
    C6half(:,j)=-(force(z,e/2)-force(z,-e/2))/h(j);
    Jq6(:,j)=(force(e,z)-force(-e,z))/(2*h(j));
    Jq6half(:,j)=(force(e/2,z)-force(-e/2,z))/h(j);
end
report.topHeight_m=top(3);
report.C4_code=rotor.dampingMatrix4;
report.C4_forceDerivative=C6(idx,idx);
report.C4_smallAngle_vs_forceDerivativeError=norm(Csmall-C6(idx,idx),'fro')/norm(C6(idx,idx),'fro');
report.C4_relativeFrobeniusError=norm(rotor.dampingMatrix4-C6(idx,idx),'fro')/norm(C6(idx,idx),'fro');
report.C4_blockRelativeErrors=zeros(2);
for i=1:2
 for j=1:2
  ii=(i-1)*2+(1:2); jj=(j-1)*2+(1:2);
  report.C4_blockRelativeErrors(i,j)=norm(rotor.dampingMatrix4(ii,jj)-C6(idx(ii),idx(jj)),'fro')/norm(C6(idx(ii),idx(jj)),'fro');
 end
end
report.velocityFD_stepHalvingError=norm(C6-C6half,'fro')/norm(C6half,'fro');
report.positionFD_stepHalvingError=norm(Jq6-Jq6half,'fro')/max(1,norm(Jq6half,'fro'));
report.meanLoad6_N_Nm=Q0;
report.loadPositionDerivative6=Jq6;
report.C6_forceDerivative=C6;
report.rotation001rad_loadChange4=zeros(4,2);
for j=1:2
 q=z; q(j+3)=0.001;
 dq=force(q,z)-Q0;
 report.rotation001rad_loadChange4(:,j)=dq(idx);
end
% Compare finite-difference full-six-channel feedback after structural projection.
model=assemble_tower_rna_structure(nodes,elems,conn,rna.massMatrix6);
free=(7:model.numDofTotal).';
modal=solve_tower_rna_undamped_modes(model.K_global(free,free),model.M_global(free,free),16);
Phi=modal.modeShapeUndamped(:,1:6);
H=model.towerTopSelectionMatrix(:,free)*Phi;
report.freq6_Hz=modal.naturalFreqUndamped_Hz(1:6);
report.Cmodal4=H(idx,:).'*rotor.dampingMatrix4*H(idx,:);
report.Cmodal6=H.'*C6*H;
report.KaeroModal=-H.'*Jq6*H;
report.KstructModal=Phi.'*model.K_global(free,free)*Phi;
report.modalDiagonalAeroStiffnessRatio=diag(report.KaeroModal)./diag(report.KstructModal);
report.aoaRange_deg=zeros(3,2);
for b=1:3
 a=rad2deg(rotor.bladeAeroOut(b).sectionAoA_rad);
 report.aoaRange_deg(b,:)=[min(a),max(a)];
end
disp(jsonencode(report));
end

function Q = independent_force(q,v,rotor,blade,op,db,sections)
% q=[global translation; rotation vector], v=[translation velocity; angular velocity].
% Uniform wind makes pure translation irrelevant; moments are about moving tower top.
a=q(4:6); S=[0 -a(3) a(2);a(3) 0 -a(1);-a(2) a(1) 0];
U=expm(S); Q=zeros(6,1);
for b=1:3
 g=rotor.bladeAeroOut(b); gp=g.gaussPointInfo;
 T=U*g.rotationReferenceToGlobal(:,1:2);
 r=U*gp.rTowerTopToSection_global;
 n=size(r,2);
 V=T.'*(op.windVelocity_global-v(1:3)-cross(repmat(v(4:6),1,n),r,1));
 i=gp.segmentIndex;
 fraction=(gp.distanceAlongBlade-blade.distanceAlongBlade(i))./ ...
    (blade.distanceAlongBlade(i+1)-blade.distanceAlongBlade(i));
 f1=section_forces(V,g.bladePitch_rad,i,blade,op,db,sections);
 f2=section_forces(V,g.bladePitch_rad,i+1,blade,op,db,sections);
 fg=T*(f1.*(1-fraction.')+f2.*fraction.');
 Q=Q+[fg;cross(r,fg,1)]*gp.integrationWeight_m;
end
end

function f=section_forces(V,pitch,i,blade,op,db,sections)
alpha=atan2(V(1,:).',V(2,:).')-pitch-blade.aeroTwist_rad(i);
alpha=mod(alpha+pi,2*pi)-pi;
[cl,cd]=query_aerofoil_section(db,sections,blade.foilSection(i),blade.thickness_pct(i),alpha);
W=sqrt(sum(V.^2,1));
f=0.5*op.airDensity*(blade.chord(i).'.*W).* ...
 [cd.'.*V(1,:)+cl.'.*V(2,:);cd.'.*V(2,:)-cl.'.*V(1,:)];
end
