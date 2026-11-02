function dxdt = pendulum_double_angle_inertiamat_rhs(x,param)

% Unpack parameters
m1 = param(1); L1 = param(2);
m2 = param(3); L2 = param(4);
g = param(5);

% Unpack state variables
th1 = x(1);
th2 = x(2);
thd1 = x(3);
thd2 = x(4);

% Compute helper variables
s = sin(th1 - th2);
c = cos(th1 - th2);
m = m1 + m2;
a = m2*L1*L2;

% Construct mass inertia matrix
M = [m*L1^2, a*c;
     a*c, m2*L2^2];

% This thing
f = [-a*s*thd2^2 - m*g*L1*sin(th1);
      a*s*thd1^2 - m2*g*L2*sin(th2)];

% Get acceleration
d2dt = M\f;

% Reshape to ODE
dxdt = [thd1;thd2;d2dt];


end
