function db = build_aerofoil_database(xfile)
% 分段线性极曲线与一致斜率；折点采用左右斜率平均（广义切线）。
% 完整360度极曲线按周期查询；非完整极曲线范围外保持端点，导数为零。
names = string(sheetnames(xfile)); foil = struct([]);
for k=1:numel(names)
    a=readmatrix(xfile,'Sheet',names(k));
    if size(a,2)<4, error('aero:PolarColumns','翼型%s需要四列。',names(k)); end
    a=a(:,1:4); a=a(~all(isnan(a),2),:);
    if any(~isfinite(a(:))), error('aero:PolarFinite','翼型%s存在缺失或非有限数据。',names(k)); end
    [alpha,ix]=sort(deg2rad(a(:,1))); values=a(ix,2:4);
    if numel(alpha)<2 || any(diff(alpha)<=0)
        error('aero:PolarGrid','翼型%s攻角需要至少两个互异点。',names(k));
    end
    token=regexp(char(names(k)),'\d+\.?\d*','match','once');
    if isempty(token), error('aero:FoilName','翼型名%s无法解析厚弦比。',names(k)); end
    periodic=abs(alpha(end)-alpha(1)-2*pi)<1e-10;
    if periodic && max(abs(values(1,:)-values(end,:)))>1e-10
        error('aero:PeriodicEndpoints','翼型%s的360度端点系数不连续。',names(k));
    end
    foil(k).name=names(k); foil(k).sheetName=names(k);
    foil(k).tc_percent=str2double(token);
    foil(k).alpha_rad=alpha; foil(k).alpha_deg=rad2deg(alpha);
    foil(k).Cl=values(:,1); foil(k).Cd=values(:,2); foil(k).Cm=values(:,3);
    foil(k).periodic=periodic;
    foil(k).query=@(q) linear_polar_query(alpha,values,q,periodic);
end
db.xfile=xfile; db.foil=foil; db.nameList=names(:);
db.tcPercentList=[foil.tc_percent].'; db.query=@(name,q) query_name(foil,name,q);
db.derivativeMethod='piecewise-linear exact slopes; mean at knots';
end
function varargout=query_name(foil,name,q)
i=find(string({foil.name})==string(name),1);
if isempty(i), error('aero:MissingFoil','未找到翼型%s。',string(name)); end
[varargout{1:nargout}]=foil(i).query(q);
end
