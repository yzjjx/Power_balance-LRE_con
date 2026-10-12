% 计算 N 轴机器人的 Rayleigh 损耗函数标量 FR。
% 运行顺序：RNEA_SYM; generate_FR;
% 粘性摩擦：FR = 1/2 * sum(Fv_i * dq_i^2)。
% 对速度求偏导得到需克服的摩擦力矩；实际摩擦作用力矩取负号。

assert(exist('N', 'var') && exist('mat_dq', 'var'), ...
    '请先运行 RNEA_SYM.m，生成 N 和 mat_dq。');
assert(isa(mat_dq, 'sym') && isequal(size(mat_dq), [N, 1]), ...
    'mat_dq 必须是 N×1 关节速度符号列向量。');

%% 1. 粘性摩擦系数：沿用工作区参数，否则生成符号参数
if ~exist('mat_Fv_value', 'var')
    mat_Fv_value = sym('Fv_', [N, 1], 'real');
end
% 数值示例（N=2 时）：取消下一行注释并填写实际非负系数。
% mat_Fv_value = [0.1; 0.2];
assert((isa(mat_Fv_value, 'sym') || isnumeric(mat_Fv_value)) && ...
    isequal(size(mat_Fv_value), [N, 1]), ...
    'mat_Fv_value 必须是 N×1 摩擦系数列向量。');

%% 2. 粘性摩擦损耗函数
FR = simplify(sym(1)/2 * sum(mat_Fv_value .* mat_dq.^2));

%% 3. 库伦摩擦扩展：取消以下三行注释即可启用
% assume(mat_dq, 'real'); % abs 的速度偏导按实数速度计算
% mat_Fc_value = sym('Fc_', [N, 1], 'real'); % 非负系数，数值设置方式同上
% FR = simplify(FR + sum(mat_Fc_value .* abs(mat_dq)));
% 此项是非光滑耗散势：非零速度处的偏导为 Fc_i*sign(dq_i)。
% 零速度处不可微，静摩擦需单独建模。

%% 4. 速度偏导给出摩擦力矩，并计算耗散功率
tau_FR = simplify(jacobian(FR, mat_dq).');
P_loss = simplify(mat_dq.' * tau_FR);
% 仅含粘性摩擦时 P_loss = 2*FR；启用库伦摩擦后按上式计算。

disp('FR =');
disp(FR);
