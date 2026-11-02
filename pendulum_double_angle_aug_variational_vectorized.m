function dx_augdt = pendulum_double_angle_aug_variational_vectorized(x_aug, param)
% Unpack states and variations
x = x_aug(1:4,:);
delta_x = x_aug(5:8,:);

% Propagate states and varaitions
dxdt = pendulum_double_angle_minTrig_rhs(x, param);
ddelta_xdt = pendulum_double_angle_minTrig_Jdx_vectorized(x, delta_x, param); % Propagate difference trajectory: δx

% Combine to a singel augmented RHS
dx_augdt = [dxdt;ddelta_xdt];

end
