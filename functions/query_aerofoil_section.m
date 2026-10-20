function [Cl,Cd,Cm,dCl,dCd,dCm]=query_aerofoil_section(db,sets,section,tc,alpha)
% 截面编号内按相邻厚弦比线性插值；按翼型分组查询，保持系数/导数一致。
counts=[numel(section),numel(tc),numel(alpha)]; n=max(counts);
if n==0 || any(counts~=1 & counts~=n), error('aero:QuerySize','查询尺寸不一致。'); end
if numel(alpha)>1, sz=size(alpha); elseif numel(tc)>1, sz=size(tc); else, sz=size(section); end
section=expand(section,n); tc=expand(tc,n); alpha=expand(alpha,n);
if any(~isfinite([section;tc;alpha])), error('aero:QueryFinite','查询输入必须有限。'); end
lo=strings(n,1); hi=lo; w=zeros(n,1);
for id=unique(section).'
    si=find([sets.sectionNumber]==id,1);
    if isempty(si), error('aero:MissingSection','未找到截面编号%g。',id); end
    ids=find(section==id); x=sets(si).tcPercentList; names=sets(si).foilNameList;
    for j=ids.'
        if numel(x)==1 || tc(j)<=x(1), l=1; u=1;
        elseif tc(j)>=x(end), l=numel(x); u=l;
        else, u=find(x>=tc(j),1); l=u-1; w(j)=(tc(j)-x(l))/(x(u)-x(l)); end
        lo(j)=names(l); hi(j)=names(u);
    end
end
needSlope=nargout>3; nv=3+3*needSlope; result=zeros(n,nv);
for name=unique([lo;hi]).'
    weight=(lo==name).*(1-w)+(hi==name).*w;
    ids=find(weight~=0); out=cell(1,nv);
    [out{:}]=db.query(name,alpha(ids));
    result(ids,:)=result(ids,:)+weight(ids).*cat(2,out{:});
end
Cl=reshape(result(:,1),sz); Cd=reshape(result(:,2),sz); Cm=reshape(result(:,3),sz);
if needSlope
 dCl=reshape(result(:,4),sz); dCd=reshape(result(:,5),sz); dCm=reshape(result(:,6),sz);
end
end
function x=expand(x,n)
x=x(:); if numel(x)==1, x=repmat(x,n,1); end
end
