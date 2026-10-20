function Q=evaluate_parked_aero_force(q,v,aero,db,sets)
% 验证探针：更新刚体姿态、相对来流、力方向和力臂，独立计算总载荷。
% q旋转采用旋转向量，v旋转分量为全局角速度；返回全局物理力/力矩。
% 仅用于零状态切线和小转角适用性诊断，不是完整非线性结构求解器。
q=q(:); v=v(:);
if numel(q)~=6 || numel(v)~=6 || any(~isfinite([q;v]))
 error('aero:ProbeInput','q,v必须为有限六元素向量。');
end
a=q(4:6); U=expm([0,-a(3),a(2);a(3),0,-a(1);-a(2),a(1),0]);
blade=aero.blade; op=aero.operatingPoint; Q=zeros(6,1);
for b=1:numel(aero.rotor4.bladeAeroOut)
 g=aero.rotor4.bladeAeroOut(b); gp=g.gaussPointInfo;
 T=U*g.rotationReferenceToGlobal(:,1:2); r=U*gp.rTowerTopToSection_global;
 V=T.'*(op.windVelocity_global(:)-v(1:3)-cross(repmat(v(4:6),1,size(r,2)),r,1));
 i=gp.segmentIndex;
 t=(gp.distanceAlongBlade-blade.distanceAlongBlade(i))./(blade.distanceAlongBlade(i+1)-blade.distanceAlongBlade(i));
 f1=section_force(V,i,g.bladePitch_rad,blade,op,db,sets);
 f2=section_force(V,i+1,g.bladePitch_rad,blade,op,db,sets);
 f=T*(f1.*(1-t.')+f2.*t.');
 Q=Q+[f;cross(r,f,1)]*gp.integrationWeight_m;
end
end
function f=section_force(V,i,pitch,blade,op,db,sets)
alpha=mod(atan2(V(1,:).',V(2,:).')-pitch-blade.aeroTwist_rad(i)+pi,2*pi)-pi;
[cl,cd]=query_aerofoil_section(db,sets,blade.foilSection(i),blade.thickness_pct(i),alpha);
f=0.5*op.airDensity*(blade.chord(i).'.*sqrt(sum(V.^2,1))).* ...
 [cd.'.*V(1,:)+cl.'.*V(2,:);cd.'.*V(2,:)-cl.'.*V(1,:)];
end
