function aero=build_parked_aero_linearization(blade,rna,top,op,db,sets)
% 未变形参考状态：Q(q,v)=Q0+Jq*q-C*v。全局力矩关于随动塔顶。
% q,v顺序[uX,uY,uZ,thetaX,thetaY,thetaZ]，不对Jq或C强制对称。
% 采用小转角增量/全局物理力矩。非线性探针仅用于验证一阶切线。
old=build_parked_rotor_aero_damping_4dof(blade,rna,top,op,db,sets);
nb=numel(old.bladeAeroOut); Qby=zeros(6,nb); Cby=zeros(6,6,nb); Jby=Cby;
for b=1:nb
    g=old.bladeAeroOut(b); gp=g.gaussPointInfo;
    R=g.rotationReferenceToGlobal; T=R(:,1:2); V=g.windVelocity_reference;
    W=hypot(V(1),V(2));
    fs=0.5*op.airDensity*(blade.chord(:).'*W).* ...
        [g.sectionDragCoeff.'*V(1)+g.sectionLiftCoeff.'*V(2); ...
         g.sectionDragCoeff.'*V(2)-g.sectionLiftCoeff.'*V(1)];
    for k=1:numel(gp.integrationWeight_m)
        i=gp.segmentIndex(k); r=gp.rTowerTopToSection_global(:,k);
        t=(gp.distanceAlongBlade(k)-blade.distanceAlongBlade(i))/ ...
          (blade.distanceAlongBlade(i+1)-blade.distanceAlongBlade(i));
        f=T*((1-t)*fs(:,i)+t*fs(:,i+1));
        c2=gp.sectionDampingMatrix2(:,:,k);
        B=T.'*[eye(3),-skew3(r)]; weight=gp.integrationWeight_m(k);
        Qby(:,b)=Qby(:,b)+[f;cross(r,f)]*weight;
        Cby(:,:,b)=Cby(:,:,b)+B.'*c2*B*weight;
        % d(T*f_local)/dtheta = -S(f)+T*C2*T'*S(U)。
        % d(r cross f)/dtheta = S(f)*S(r)+S(r)*df/dtheta。
        df=-skew3(f)+T*c2*T.'*skew3(op.windVelocity_global(:));
        dm=skew3(f)*skew3(r)+skew3(r)*df;
        Jby(:,4:6,b)=Jby(:,4:6,b)+[df;dm]*weight;
    end
end
aero.Q0=sum(Qby,2); aero.C6=sum(Cby,3); aero.Jq6=sum(Jby,3);
aero.Ka6=-aero.Jq6; aero.Q0_byBlade=Qby; aero.C6_byBlade=Cby; aero.Jq6_byBlade=Jby;
aero.rotor4=old; aero.blade=blade; aero.operatingPoint=op;
aero.referenceDescription='undeformed geometry; uniform wind; infinitesimal global rotations';
end
function S=skew3(r)
S=[0,-r(3),r(2);r(3),0,-r(1);-r(2),r(1),0];
end
