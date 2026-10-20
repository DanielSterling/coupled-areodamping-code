function [Cl,Cd,Cm,dCl,dCd,dCm]=linear_polar_query(alpha,values,q,periodic)
% values每行为[Cl,Cd,Cm]，alpha单位rad。折点通常不可微，返回平均切线。
sz=size(q); x=q(:);
if any(~isfinite(x)), error('aero:QueryFinite','查询攻角必须有限。'); end
if periodic, x=mod(x-alpha(1),2*pi)+alpha(1); end
outside=x<alpha(1) | x>alpha(end); xc=min(max(x,alpha(1)),alpha(end));
bin=discretize(xc,alpha); bin=min(bin,numel(alpha)-1);
slopes=diff(values)./diff(alpha);
v=values(bin,:)+(xc-alpha(bin)).*slopes(bin,:);
Cl=reshape(v(:,1),sz); Cd=reshape(v(:,2),sz); Cm=reshape(v(:,3),sz);
if nargout<=3, return; end
s=slopes(bin,:); s(outside,:)=0; tol=32*eps(max(1,max(abs(alpha))));
for k=1:numel(alpha)
    hit=abs(x-alpha(k))<=tol;
    if ~any(hit), continue; end
    if k==1 || k==numel(alpha)
        if periodic, sk=(slopes(1,:)+slopes(end,:))/2;
        elseif k==1, sk=slopes(1,:)/2;
        else, sk=slopes(end,:)/2; end
    else, sk=(slopes(k-1,:)+slopes(k,:))/2; end
    s(hit,:)=repmat(sk,sum(hit),1);
end
dCl=reshape(s(:,1),sz); dCd=reshape(s(:,2),sz); dCm=reshape(s(:,3),sz);
end
