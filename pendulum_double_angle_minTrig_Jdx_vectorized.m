function ddelta_xdt = pendulum_double_angle_minTrig_Jdx_vectorized(x, delta_x, param)
% Computes the right-hand side of the variational part of the augmented ODE
% It solves this part: dδx/dt = J*δx

% State
alpha     = x(1,:);
beta      = x(2,:);
dalpha_dt = x(3,:);
dbeta_dt  = x(4,:);

% Variation
dxa  = delta_x(1,:);
dxb  = delta_x(2,:);
dxad = delta_x(3,:);
dxbd = delta_x(4,:);

% Parameters
m1 = param(1,:);
L1 = param(2,:);
m2 = param(3,:);
L2 = param(4,:);
g  = param(5,:);

% These 6 trig functions are all we need, everything else is algebra
cos_a = cos(alpha);
sin_a = sin(alpha);
cos_b = cos(beta);
sin_b = sin(beta);
cos_a_b = cos(alpha - beta);
sin_a_b = sin(alpha - beta);

% Using some identities to avoid further calls to trig functions
cos_a_b_sq = cos_a_b.^2;
sin_a_b_sq = 1 - cos_a_b_sq;                  % 1 - cos^2(α−β)
cos_2a_2b = 2.*cos_a_b_sq - 1;                % cos(2(α−β))
sin_2a_2b = 2.*sin_a_b.*cos_a_b;              % sin(2(α−β))
cos_2a_b = cos_a.*cos_a_b - sin_a.*sin_a_b;   % cos(2α−β))
sin_2a_b = sin_a.*cos_a_b + cos_a.*sin_a_b;   % sin(2α−β))

% Further precomputations
t8 = dalpha_dt.^2;
t9 = dbeta_dt.^2;
t10 = 1.0./L1;
t11 = 1.0./L2;
t16 = m2./2.0;
t14 = g.*m1.*sin_a;
t15 = g.*m2.*sin_a;
t20 = (g.*m1.*sin_b)./2.0;
t25 = g.*m2.*sin_b.*(-1.0./2.0);
t24 = -t20;
t33 = g.*m2.*sin_b.*cos_a_b;
t34 = g.*m2.*sin_b.*sin_a_b;
t36 = L1.*m1.*t8.*cos_a_b;
t37 = L1.*m2.*t8.*cos_a_b;
t38 = L2.*m2.*t9.*cos_a_b;
t39 = L1.*m1.*t8.*sin_a_b;
t40 = L1.*m2.*t8.*sin_a_b;
t41 = L2.*m2.*t9.*sin_a_b;
t48 = (g.*m1.*sin_2a_b)./2.0;
t49 = g.*t16.*sin_2a_b;
t32 = m2.*cos_a_b_sq;
t43 = -t33;
t44 = m2.*cos_2a_2b.*(-1.0./2.0);
t45 = L2.*m2.*t9.*cos_2a_2b;
t47 = L1.*m2.*t8.*sin_a_b_sq;
t52 = L2.*t9.*t16.*sin_2a_2b;
t55 = sin_a_b.*t37;
t35 = -t32;
t46 = L1.*t8.*t32;
t51 = -t47;
t56 = m1+t16+t44;
t59 = t14+t15+t41+t43+t55;
t60 = t24+t25+t39+t40+t48+t49+t52;
t50 = m1+m2+t35;
t57 = 1.0./t56;
t53 = 1.0./t50;
t58 = t57.^2;
t54 = t53.^2;
t62 = m2.*t11.*sin_2a_2b.*t58.*t60;
t61 = m2.*t10.*cos_a_b.*sin_a_b.*t54.*t59.*2.0;

% Compute each element of the Jacobian individually since they wont be scalars but 1 x Np vectors where each element corresponds one Jacobian value for a given pendulum. Like: [J_ij1, J_ij2, ... J_ijNp]
% These are the elemnets in the first column of the Jacobian: J_i1 = df_i/dx1 
% J_11 = 0;
% J_21 = 0;
J_31 = t61-t10.*t53.*(t34+t38+t46+t51+g.*m1.*cos_a+g.*m2.*cos_a);
J_41 = -t62+t11.*t57.*(t36+t37+t45+g.*m1.*cos_2a_b+g.*m2.*cos_2a_b);
% Second column of Jacobian: J_i2 = df_i/dx2 
% J_12 = 0;
% J_22 = 0;
J_32 = -t61+t10.*t53.*(t34+t38+t46+t51+g.*m2.*cos_b.*cos_a_b);
J_42 = t62-t11.*t57.*(t36+t37+t45+(g.*m1.*cos_b)./2.0+(g.*m1.*cos_2a_b)./2.0+g.*cos_b.*t16+g.*t16.*cos_2a_b);
% Third column of Jacobian: J_i3 = df_i/dx3
% J_13 = 1;
% J_23 = 0;
J_33 = dalpha_dt.*m2.*cos_a_b.*sin_a_b.*t53.*-2.0;
J_43 = t11.*t57.*(L1.*dalpha_dt.*m1.*sin_a_b.*2.0+L1.*dalpha_dt.*m2.*sin_a_b.*2.0);
% Fourth column of Jacobian: J_i4 = df_i/dx4
% J_14 = 0;
% J_24 = 1;
J_34 = L2.*dbeta_dt.*m2.*t10.*sin_a_b.*t53.*-2.0;
J_44 = dbeta_dt.*m2.*sin_2a_2b.*t57;

% J = [0, 0, 1, 0;
%      0, 0, 0, 1;
%      J_31, J_32, J_33, J_34;
%      J_41, J_42, J_43, J_44 ];

% Compute dδx/dt = J*δx suitable for vectorized inputs. So dxa is a 1 x Np vector where Np is the number of simultaneously simulated double pendulums. 
ddelta_xdt = zeros(4, size(delta_x,2), 'like', delta_x);
ddelta_xdt(1,:) = dxad;
ddelta_xdt(2,:) = dxbd;
ddelta_xdt(3,:) = J_31.*dxa + J_32.*dxb + J_33.*dxad + J_34.*dxbd;
ddelta_xdt(4,:) = J_41.*dxa + J_42.*dxb + J_43.*dxad + J_44.*dxbd;


end









