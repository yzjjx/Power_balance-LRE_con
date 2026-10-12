# 功率回归平衡方程推导
该推导依赖MATALAB二连杆相关符号动力学推导代码：https://github.com/yzjjx/min_parameter  

在此基础上，在tau中加入摩擦力
## 1 粘性摩擦与库伦摩擦
黏性摩擦力矩与关节转速近似成正比：

$$
\tau_v=f_v\dot q
$$

其中 $f_v\ge0$ 为黏性摩擦系数。对于旋转关节，其单位为 N·m·s/rad。

*************************************

理想库仑摩擦力矩的幅值不随非零转速变化，仅与转动方向有关：

$$
\tau_c=f_c\operatorname{sgn}(\dot q)
$$

其中 $f_c\ge0$ 为库仑摩擦幅值（N·m），$\operatorname{sgn}$ 为符号函数：

$$
\operatorname{sgn}(v)=\begin{cases}
+1,&v>0,\\
0,&v=0,\\
-1,&v<0.
\end{cases}
$$

因此在采用 $\operatorname{sgn}(0)=0$ 的简单数值模型时：

$$
\tau_c=\begin{cases}
+f_c,&\dot q>0,\\
0,&\dot q=0,\\
-f_c,&\dot q<0.
\end{cases}
$$

**注意静止状态**：在理想库仑摩擦的静力学解释中，$\dot q=0$ 时可能存在使关节保持静止的静摩擦力矩，其取值通常是一个区间，而非必然为零。例如可用集合值形式近似表达：

$$
\tau_{\mathrm{static}}\in[-f_s,f_s]\qquad(\dot q=0)
$$

其中 $f_s$ 是最大静摩擦力矩。$\operatorname{sgn}(0)=0$ 是编程上的一种单值选择，不能描述完整静摩擦与黏滑行为。

*************************
最终加入的MATLAB文件为：MATLAB\add_friction.m，仅加入粘性摩擦，取消注释后加入库伦摩擦

## 2 重力势能提取

最终加入的MATLAB文件为：MATLAB\generate_U.m，用来输出势能函数

## 4 加入一阶滤波器
