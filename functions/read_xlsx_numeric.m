function [X, headers] = read_xlsx_numeric(xlsxFile, sheetName)
% [X, headers] = read_xlsx_numeric(xlsxFile, sheetName)
% 读取 xlsx，返回：
%   headers : 1×n 的列名 string
%   X       : 仅保留“整行均可转为数值”的纯数值矩阵
%
% 规则：
%   1) 数值单元格：直接保留
%   2) 字符串形式的数字：尝试转成数值
%   3) 只要某一行存在无法转成数值的单元格，则整行删除

    if nargin < 2
        sheetName = 1;
    end

    T = readtable(xlsxFile, 'Sheet', sheetName, 'PreserveVariableNames', true);

    headers = string(T.Properties.VariableNames);

    nRow = height(T);
    nCol = width(T);

    % 先构造一个同尺寸的数值矩阵
    X_all = nan(nRow, nCol);
    validRow = true(nRow, 1);   % 标记这一行是否全都可转成数值

    for j = 1:nCol
        col = T{:, j};

        if isnumeric(col)
            xj = double(col);

        elseif iscell(col) || isstring(col) || ischar(col)
            s = string(col);
            s = strtrim(s);
            xj = str2double(s);

        else
            % 其他类型（datetime、logical、categorical 等）整列都视为不可用
            xj = nan(nRow, 1);
        end

        X_all(:, j) = xj;

        % 这一列里如果某行转不成数值，则该行删除
        validRow = validRow & ~isnan(xj);
    end

    % 只保留整行均为数值的行
    X = X_all(validRow, :);
end