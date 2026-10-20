function impact=read_impact_load(path,startTime)
% 严格读取NUMPOINTS和两列数据；零端点保证文件范围外零载荷无伪斜坡。
fid=fopen(path,'r');
if fid<0, error('load:MissingFile','无法读取Impact文件：%s',path); end
cleanup=onCleanup(@() fclose(fid));
header=strtrim(fgetl(fid));
token=regexp(header,'^NUMPOINTS\s+(\d+)\s*$','tokens','once');
if isempty(token), error('load:Header','Impact首行必须为NUMPOINTS N；请检查是否被加密。'); end
lines=splitlines(string(fread(fid,Inf,'*char').'));
data=zeros(0,2);
for k=1:numel(lines)
 line=strtrim(lines(k)); if strlength(line)==0, continue; end
 items=split(line); items=items(strlength(items)>0);
 if numel(items)~=2, error('load:Columns','Impact每个数据行必须为两列。'); end
 data(end+1,:)=str2double(items).'; %#ok<AGROW>
end
if size(data,1)~=str2double(token{1}) || size(data,1)<2
 error('load:Count','Impact实际数据点数与NUMPOINTS不符或不足两点。');
end
if any(~isfinite(data(:))) || any(diff(data(:,1))<=0)
 error('load:Time','Impact数据必须有限，时间必须严格递增。');
end
tol=1e-12*max(1,max(abs(data(:,2))));
if any(abs(data([1,end],2))>tol)
 error('load:Endpoints','Impact首尾力必须为零；非零端点需另行定义不连续加载规则。');
end
data([1,end],2)=0;
impact.timeLocal=data(:,1); impact.timeGlobal=startTime+data(:,1); impact.force=data(:,2);
impact.path=path; impact.startTime=startTime;
end
