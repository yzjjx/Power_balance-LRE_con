% 使用惯性矩阵 M 构建 N 轴机器人的动能标量 T。
% 运行顺序：RNEA_SYM; Mat_out; generate_T;
% Mat_out 已通过 M = jacobian(tau, mat_ddq) 提取惯性矩阵。
% 动能公式：T(q,dq) = 1/2 * dq.' * M(q) * dq。

assert(exist('M', 'var') && exist('mat_dq', 'var') && exist('N', 'var'), ...
    '请先运行 RNEA_SYM.m 和 Mat_out.m，生成 M、mat_dq 和 N。');
assert((isa(M, 'sym') || isnumeric(M)) && isequal(size(M), [N, N]) && ...
    isa(mat_dq, 'sym') && isequal(size(mat_dq), [N, 1]), ...
    'M 必须是 N×N 惯性矩阵，mat_dq 必须是 N×1 关节速度符号列向量。');

%% 计算动能：结果为标量，速度为零时动能为零
T = simplify(sym(1)/2 * mat_dq.' * M * mat_dq);

disp('T =');
disp(T);
