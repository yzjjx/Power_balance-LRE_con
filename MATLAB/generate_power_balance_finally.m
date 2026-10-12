% 构建一阶滤波后的功率回归：y = Omega.'*theta，适用于 N 轴机器人。
% 运行顺序：RNEA_SYM; Mat_out; generate_T; generate_U; generate_FR;
%           generate_power_balance_finally;
% 本脚本自动调用 generate_theta_power_bal，以相同的 SVD + 两次 QR 提取 theta。
% H(p)=1/(p+lambda)：能量基函数使用 H(p)[d(phi_E)/dt]，摩擦功率使用 H(p)[phi_D]。
% 注意：论文式(68)将 pH(p)作用于未积分的摩擦功率，与当前粘性摩擦模型不一致。
% 此处保留能量平衡正确的形式；不对摩擦功率额外求导。

%% 1. 提取参数，并将相对机械能、摩擦功率写成同一个 theta 的线性组合
generate_theta_power_bal;
power_energy_origin = subs(T + U, [mat_q; mat_dq], sym(zeros(2*N, 1)));
power_energy_relative = simplify(T + U - power_energy_origin);
% 移除与运动无关的能量常数，避免它进入参数辨识；不修改原 T、U。
power_basis_models = [power_energy_relative; loss_power];
power_basis_full = sym(zeros(2, numel(theta_power_full)));
for model_idx = 1:2
    [coefficients, terms] = coeffs(expand(power_basis_models(model_idx)), power_parameter_symbols);
    for term_idx = 1:numel(terms)
        [found, column_idx] = ismember(char(terms(term_idx)), power_term_names);
        if ~found && ~isequal(coefficients(term_idx), sym(0))
            error('generate_power_balance_finally:MissingParameterTerm', ...
                '能量或摩擦功率包含未提取的参数项，请检查前置模型。');
        elseif found
            power_basis_full(model_idx, column_idx) = ...
                power_basis_full(model_idx, column_idx) + coefficients(term_idx);
        end
    end
end
power_basis = simplify(power_basis_full(:, power_base_columns));
phi_E = power_basis(1, :).'; % 相对机械能 = phi_E.'*theta
phi_D = power_basis(2, :).'; % 摩擦耗散功率 = phi_D.'*theta
power_basis_residual = simplify(expand(power_basis_full - power_basis*power_parameter_map));
assert(all(isAlways(power_basis_residual == 0), 'all') && ...
    all(isAlways(simplify(expand(power_basis_models - power_basis*theta)) == 0), 'all'), ...
    '无法使用同一个 theta 表示机械能与摩擦功率。');

%% 2. 一阶滤波器的连续时间状态方程，默认 lambda=1（论文实验设置）
lambda = 1; % 可修改，必须为正数
assert(isnumeric(lambda) && isscalar(lambda) && isfinite(lambda) && lambda > 0, ...
    'lambda 必须是正的有限数值标量。');
power_w = numel(theta);
z_E = sym('z_E_', [power_w, 1], 'real');
z_D = sym('z_D_', [power_w, 1], 'real');
y = sym('y_power', 'real');
tau_power_input = sym('tau_power_', [N, 1], 'real');
% tau_power_input 表示实测总驱动力矩（包含克服摩擦所需的力矩）。
% 仅 RNEA_SYM 生成的无摩擦 tau 不能作为含摩擦模型的总驱动力矩。
power_input = mat_dq.' * tau_power_input;
dz_E = -lambda*z_E + phi_E;
dz_D = -lambda*z_D + phi_D;
dy = -lambda*y + power_input;
Omega = simplify(phi_E - lambda*z_E + z_D);
power_regression_equation = y == Omega.'*theta;
power_filter_state = [z_E; z_D; y];
power_filter_rhs_symbolic = [dz_E; dz_D; dy];

%% 3. 初值：保证 y(0)=0、Omega(0)=0，兼容非零初始姿态和速度
power_q0_symbols = sym('q0_power_', [N, 1], 'real');
power_dq0_symbols = sym('dq0_power_', [N, 1], 'real');
phi_E0 = subs(phi_E, [mat_q; mat_dq], [power_q0_symbols; power_dq0_symbols]);
power_filter_initial_symbolic = [phi_E0/lambda; sym(zeros(power_w + 1, 1))];
% z_E(0)=phi_E(0)/lambda，使 phi_E-lambda*z_E 严格实现 H(p)[d(phi_E)/dt]。
% z_D(0)=0，y(0)=0；无需计算或测量关节加速度。

%% 4. 验证：r=y-Omega.'*theta 满足 dr/dt+lambda*r=输入功率-模型功率
dot_Omega = jacobian(phi_E, mat_q)*mat_dq + ...
            jacobian(phi_E, mat_dq)*mat_ddq - lambda*dz_E + dz_D;
power_filter_identity_residual = simplify(expand(dy - dot_Omega.'*theta + ...
    lambda*(y - Omega.'*theta) - (power_input - power_model)));
power_filter_initial_residual = simplify(subs(y - Omega.'*theta, ...
    [power_filter_state; mat_q; mat_dq], ...
    [power_filter_initial_symbolic; power_q0_symbols; power_dq0_symbols]));
assert(isAlways(power_filter_identity_residual == 0) && ...
    isAlways(power_filter_initial_residual == 0), '滤波回归的状态方程或初值验证失败。');

%% 5. 数值接口：已知参数取值与 generate_theta_power_bal 的采样设置一致
% 只替换 g 等已知量，不代入未知 theta；接口中不使用 ddq。
power_phi_E_numeric = simplify(subs(phi_E, power_known_symbols, power_known_values));
power_phi_D_numeric = simplify(subs(phi_D, power_known_symbols, power_known_values));
power_remaining_symbols = setdiff(symvar([power_phi_E_numeric; power_phi_D_numeric]).', ...
    [mat_q; mat_dq], 'stable');
assert(isempty(power_remaining_symbols), ...
    '基函数仍含其他参数，请在 power_known_symbols/values 中填写已知量。');
power_phi_E_fun = matlabFunction(power_phi_E_numeric, 'Vars', {mat_q, mat_dq});
power_phi_D_fun = matlabFunction(power_phi_D_numeric, 'Vars', {mat_q, mat_dq});
% x 的排列：[z_E; z_D; y]。输入 q、dq、tau_total 均为 N×1 列向量。
power_filter_init = @(q0,dq0) [power_phi_E_fun(q0,dq0)/lambda; zeros(power_w+1,1)];
power_filter_rhs = @(x,q,dq,tau_total) [ ...
    -lambda*x(1:power_w) + power_phi_E_fun(q,dq); ...
    -lambda*x(power_w+1:2*power_w) + power_phi_D_fun(q,dq); ...
    -lambda*x(end) + dq.'*tau_total];
power_filter_output = @(x,q,dq) deal(x(end), ...
    power_phi_E_fun(q,dq) - lambda*x(1:power_w) + x(power_w+1:2*power_w));
% 用法：x0=power_filter_init(q0,dq0);
% 连续积分：ode45(@(t,x) power_filter_rhs(x,q_fun(t),dq_fun(t),tau_fun(t)),tspan,x0);
% 对积分得到的状态 x，调用 [y_value,Omega_value]=power_filter_output(x,q,dq)。
% 实际滤波结果需用采样信号或轨迹函数积分上述状态方程。

disp('最终功率回归方程 y = Omega^T*theta：');
disp(power_regression_equation);
disp('Omega ='); disp(Omega);
disp('滤波初值 ='); disp(power_filter_initial_symbolic);
disp('验证通过：初值一致，且滤波回归满足功率平衡。');
