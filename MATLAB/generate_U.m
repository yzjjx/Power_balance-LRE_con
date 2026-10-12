% 从工作区重力向量 G 积分得到重力势能标量 U，适用于 N 轴机器人。
% 运行顺序：RNEA_SYM; Mat_out; generate_U;
% 方法一：逐关节积分，再用其余偏导关系补齐未知函数。
% 不修改 G 和 mat_q；以所有关节坐标为零的构型作为势能零点。

assert(exist('G', 'var') && exist('mat_q', 'var') && exist('N', 'var'), ...
    '请先运行 RNEA_SYM.m 和 Mat_out.m，生成 G、mat_q 和 N。');
assert(isa(G, 'sym') && isa(mat_q, 'sym') && ...
    isequal(size(G), [N, 1]) && isequal(size(mat_q), [N, 1]), ...
    'G 和 mat_q 必须是 N×1 符号列向量，关节坐标应保留为符号。');

%% 1. 从最后一个关节开始积分，逐步补齐势能
U = sym(0);
for i = N:-1:1
    % G(i) 中尚未被当前 U 的偏导表示的部分，就是待补齐函数的偏导。
    U_remainder = simplify(G(i) - diff(U, mat_q(i)));
    U = simplify(U + int(U_remainder, mat_q(i)));
end

% %% 2. 确定积分常数，使 U(0,...,0) = 0
% U = simplify(U - subs(U, mat_q, sym(zeros(N, 1))));

%% 3. 反向求偏导验证：dU/dq 必须等于 G
G_from_U = simplify(jacobian(U, mat_q).');
U_residual = simplify(expand(G_from_U - G));
assert(all(isAlways(U_residual == 0), 'all'), ...
    'generate_U:GradientMismatch', ...
    '势能偏导与 G 不一致，请检查重力向量或积分结果。');

disp('U =');
disp(U);
