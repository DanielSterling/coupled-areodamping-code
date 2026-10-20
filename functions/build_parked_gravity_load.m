function gravity=build_parked_gravity_load(structure,rna,g)
%BUILD_PARKED_GRAVITY_LOAD Consistent, frozen reference gravity load.
% A constant translational acceleration with zero nodal angular acceleration
% reproduces uniform body force through M*a. RNA offset couplings yield its
% eccentric weight moment about the tower top without another axis shift.
% No gravity tangent or axial-load geometric stiffness is introduced here.
validateattributes(g,{'numeric'},{'scalar','real','finite','nonnegative'});
if structure.dofPerNode~=6
 error('gravity:DOF','Expected six DOFs per tower node.');
end
a=zeros(structure.numDofTotal,1); a(3:6:end)=-g;
Ft=structure.M_tower_global*a;
Fr=structure.M_RNA_global*a;
P=structure.towerTopSelectionMatrix;
force=[0;0;-rna.mass*g];
expected=[force;cross(rna.rTowerTopToMassCenter(:),force)];
actual=P*Fr;
err=norm(actual-expected)/max(1,norm(expected));
if err>1e-10
 error('gravity:RNAReference','RNA gravity inconsistent with mass/COM: %.3g.',err);
end
gravity.accelerationGlobal=[0;0;-g];
gravity.accelerationField=a;
gravity.forceTowerGlobal=Ft; gravity.forceRNAGlobal=Fr;
gravity.forceGlobal=Ft+Fr;
gravity.forceRNATop=actual;
gravity.expectedRNATop=expected;
gravity.rnaReferenceRelativeError=err;
gravity.description='Uniform reference gravity; no geometric stiffness or gravity load tangent.';
end
