function dx_augdt = pendulum_double_angle_aug_variational(x_aug, param)

x = x_aug(1:4, :); % State variables [th1, th2, th1d, th2d]'
delta_x = x_aug(5:8, :); % δx

dxdt = pendulum_double_angle_minTrig_rhs(x, param);
ddelta_xdt = pendulum_double_angle_minTrig_J(x, param)*delta_x; % Propagate difference trajectory: δx

dx_augdt = [dxdt;ddelta_xdt];

end





