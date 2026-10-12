% 从 T、U、FR 提取功率平衡方程的固定参数向量 theta，适用于 N 轴机器人。
% 运行顺序：RNEA_SYM; Mat_out; generate_T; generate_U; generate_FR;
%           generate_theta_power_bal;
% 方法：随机运动状态采样 -> SVD 判秩 -> 两次 QR -> 精确符号参数组合。
% 输出：power_model = Y_power*theta（未滤波的功率平衡右端）。
% T、U、FR 和运动状态符号保持不变；不覆盖力矩回归的 Y、pi_param。

%% 1. 构建功率模型：d(T+U)/dt + dq.'*(FR 对速度的梯度)
required_power_vars = {'T','U','FR','N','mat_q','mat_dq','mat_ddq'};
for k = 1:numel(required_power_vars)
    assert(exist(required_power_vars{k}, 'var') == 1, ...
        '工作区缺少 %s，请先运行前置脚本。', required_power_vars{k});
end
assert(isa(T, 'sym') && isscalar(T) && isa(U, 'sym') && isscalar(U) && ...
    isa(FR, 'sym') && isscalar(FR), 'T、U、FR 必须是符号标量。');
assert(isa(mat_q, 'sym') && isequal(size(mat_q), [N, 1]) && ...
    isa(mat_dq, 'sym') && isequal(size(mat_dq), [N, 1]) && ...
    isa(mat_ddq, 'sym') && isequal(size(mat_ddq), [N, 1]), ...
    'q、dq、ddq 必须保留为 N×1 符号列向量。');

energy_power = T + U;
loss_power = simplify(mat_dq.' * jacobian(FR, mat_dq).');
% 粘性摩擦时 loss_power = 2*FR；此写法也支持 FR 中的库伦摩擦扩展。
power_model = simplify(jacobian(energy_power, mat_q)*mat_dq + ...
                       jacobian(energy_power, mat_dq)*mat_ddq + loss_power);
power_motion_symbols = [mat_q; mat_dq; mat_ddq];

%% 2. 设置与随机采样：只代入运动状态，保留待辨识参数
power_known_symbols = sym('g'); % g 默认已知，最终 Y_power 仍保留符号 g
power_known_values = 9.81;
% 若连杆长度已知，可改为以下设置（数值按实际模型填写）：
% power_known_symbols = [g; mat_L];
% power_known_values = [9.81; 0.45; 0.45];
power_rank_tol = 1e-10;
power_parameter_symbols = setdiff(symvar(power_model).', ...
    [power_motion_symbols; power_known_symbols], 'stable');
assert(~isempty(power_parameter_symbols), '功率模型中没有符号形式的待辨识参数。');

power_n_samples = max(20, 4*numel(power_parameter_symbols));
power_sample_stream = RandStream('mt19937ar', 'Seed', 0);
% 半角正切取有理数，使 sin(q)、cos(q) 的采样值可精确计算。
power_q_samples = 2*atan(sym(randi(power_sample_stream, [-100, 100], ...
    N, power_n_samples))/20);
power_dq_samples = sym(randi(power_sample_stream, [-6, 6], N, power_n_samples))/2;
power_ddq_samples = sym(randi(power_sample_stream, [-10, 10], N, power_n_samples))/2;
power_motion_samples = [power_q_samples; power_dq_samples; power_ddq_samples];
power_samples = sym(zeros(power_n_samples, 1));
for sample_idx = 1:power_n_samples
    power_samples(sample_idx) = simplify(subs(power_model, power_motion_symbols, ...
                                             power_motion_samples(:, sample_idx)));
end

%% 3. 从采样功率收集参数单项式，构建 W_symbolic*theta_power_full
theta_power_full = sym(zeros(0, 1));
W_power_symbolic = sym(zeros(power_n_samples, 0));
power_term_names = {};
for sample_idx = 1:power_n_samples
    [coefficients, terms] = coeffs(expand(power_samples(sample_idx)), power_parameter_symbols);
    for term_idx = 1:numel(terms)
        [found, column_idx] = ismember(char(terms(term_idx)), power_term_names);
        if ~found
            column_idx = numel(power_term_names) + 1;
            power_term_names{column_idx} = char(terms(term_idx)); %#ok<SAGROW>
            theta_power_full(column_idx, 1) = terms(term_idx);
            W_power_symbolic(:, column_idx) = sym(zeros(power_n_samples, 1));
        end
        W_power_symbolic(sample_idx, column_idx) = ...
            W_power_symbolic(sample_idx, column_idx) + coefficients(term_idx);
    end
end

%% 4. SVD 判秩，第一次 QR 选独立列，第二次 QR 选独立采样行
W_power = double(subs(W_power_symbolic, power_known_symbols, power_known_values));
power_column_scale = vecnorm(W_power, 2, 1);
power_column_scale(power_column_scale == 0) = 1;
W_power_scaled = W_power ./ power_column_scale;
power_singular_values = svd(W_power_scaled);
assert(~isempty(power_singular_values) && power_singular_values(1) > 0, ...
    '采样功率全为零，请检查输入或增加采样次数。');
power_n_base = sum(power_singular_values > power_rank_tol*power_singular_values(1));

[~, ~, power_column_order] = qr(W_power_scaled, 'vector');
power_base_columns = sort(power_column_order(1:power_n_base));
[~, ~, power_row_order] = qr(W_power_scaled(:, power_base_columns).', 'vector');
power_exact_rows = power_row_order(1:power_n_base);
power_A_exact = W_power_symbolic(power_exact_rows, :);
% 数值运算只选主元；参数关系由原符号矩阵精确求解。
power_parameter_map = simplify(power_A_exact(:, power_base_columns) \ power_A_exact);
theta = simplify(power_parameter_map * theta_power_full);

%% 5. 回到原始功率表达式，反推符号回归行向量 Y_power
Y_power_full = sym(zeros(1, numel(theta_power_full)));
[coefficients, terms] = coeffs(expand(power_model), power_parameter_symbols);
for term_idx = 1:numel(terms)
    [found, column_idx] = ismember(char(terms(term_idx)), power_term_names);
    if ~found && ~isequal(coefficients(term_idx), sym(0))
        error('generate_theta_power_bal:InsufficientSamples', ...
            '采样遗漏了功率模型中的参数项，请增加 power_n_samples。');
    elseif found
        Y_power_full(column_idx) = Y_power_full(column_idx) + coefficients(term_idx);
    end
end
Y_power_full = simplify(Y_power_full);
Y_power = Y_power_full(:, power_base_columns);

%% 6. 验证 theta 不含运动状态，且能还原任意状态下的功率模型
power_regressor_residual = simplify(expand(Y_power_full - Y_power*power_parameter_map));
power_residual = simplify(expand(power_model - Y_power*theta));
if any(has(theta, power_motion_symbols), 'all') || ...
        ~all(isAlways(power_regressor_residual == 0), 'all') || ...
        ~isAlways(power_residual == 0)
    error('generate_theta_power_bal:VerificationFailed', ...
        '符号验证失败，请增加采样次数或调整 power_rank_tol。');
end

fprintf('采样次数：%d；参数单项式：%d；独立参数组合：%d。\n', ...
    power_n_samples, numel(theta_power_full), power_n_base);
disp('功率平衡参数向量 theta =');
disp(theta);
disp('未滤波的功率回归行向量 Y_power =');
disp(Y_power);
disp('验证通过：power_model - Y_power*theta = 0。');
