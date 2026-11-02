clear
close all


%%% For a SINGLE pendulum in CARTESIAN coordinates with the pivot at the origin
syms x1(t) y1(t) m1 L1 lambda1 g t real
syms x1dd y1dd real
syms x y xd yd real


% Construct symbolic position, velocity and acceleration vectors
r = [x1(t);y1(t)]; 

x1_dot = diff(x1(t),t);
y1_dot = diff(y1(t),t);

r_dot = [x1_dot;y1_dot];

% Define kinetic and potential energy
E_kin = 1/2*m1*(x1_dot^2 + y1_dot^2);  % Same as  1/2*m1*(r_dot.'*r_dot);
E_pot = m1*g*y1(t);

% Define constraint
phi = x1(t)^2 + y1(t)^2 - L1^2;

% Lagrangian
L = E_kin - E_pot + lambda1*phi;

% % Calculate gradients
dLdr = [diff(L,x1(t));diff(L,y1(t))]; 
dLdr_dot = [diff(L,x1_dot); diff(L,y1_dot)];
ddLdr_dotdt = diff(dLdr_dot,t);

% Euler-Lagrange equation
EuLa_eq = ddLdr_dotdt - dLdr; 

% Substitute second derivatives with symbolic variables so we can solve for them
EuLa_eq_sub = subs(EuLa_eq, ...
    [diff(x1(t),t,2), diff(y1(t),t,2)], ...
    [x1dd, y1dd]);
% % Extract each line and set to zero
% eq1 = EuLa_eq_sub(1) == 0;
% eq2 = EuLa_eq_sub(2) == 0;
% Solve for the second derivatives
sol = solve(EuLa_eq_sub == 0, [x1dd, y1dd]);
xdd = sol.x1dd;
ydd = sol.y1dd;

% Group into vector
r_dot2 = [xdd;ydd];

% Use the second derivative of the pendulum link constraint
Cnstr_eq = r.'*r_dot2 + r_dot.'*r_dot;
% Solve for lambda
lambda1_sol = solve(Cnstr_eq == 0, lambda1);

% Substitute lambda back into the SECOND order ODE
d2adt = simplify(subs(r_dot2, lambda1, lambda1_sol));
% Final set of FIRST order ODE's 
dxdt = [r_dot;d2adt];


%%% Create now a .m function to use in a ODE solver after passing anonym
% First substitute the derivatives then the states. This order seems to matter. Matlab...
dxdt_sym = subs(dxdt, ...
   [ diff(x1(t),t), diff(y1(t),t)], ...
   [      xd,            yd]);

dxdt_sym = subs(dxdt_sym, ...
   [x1(t), y1(t)], ...
   [x,     y    ]);

% Group the variables
x = [x; y; xd; yd];
theta = [m1; L1; g];
% Write the function
matlabFunction(dxdt_sym, ...
    'File', 'pendulum_cartesian_rhs', ...
    'Vars', { x, theta });



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% For a DOUBLE pendulum in CARTESIAN coordinates with the pivot at the origin
syms x1(t) y1(t) x2(t) y2(t) m1 m2 L1 L2 lambda1 lambda2 g t real
syms x1dd y1dd x2dd y2dd real % used to inter-store the second derivatives 
syms X1 Y1 X1d Y1d X2 Y2 X2d Y2d real % For the final ODE's writen to .m

% Construct symbolic position, velocity and acceleration vectors
r1 = [x1(t); y1(t)]; 
r2 = [x2(t); y2(t)]; 
s = r2 - r1;

x1_dot = diff(x1(t),t);
y1_dot = diff(y1(t),t);
x2_dot = diff(x2(t),t);
y2_dot = diff(y2(t),t);

r_dot = [x1_dot; y1_dot; x2_dot; y2_dot]; 
r1_dot = [x1_dot; y1_dot];
r2_dot = [x2_dot; y2_dot];

% Define kinetic and potential energy
E_kin = 1/2*m1*(x1_dot^2 + y1_dot^2)  +  1/2*m2*(x2_dot^2 + y2_dot^2);
E_pot = m1*g*(y1(t) + L1) + m2*g*(y2(t) + L1 + L2);

% Define constraints
phi1 = x1(t)^2 + y1(t)^2 - L1^2; % Enforces constant lenght of first link
phi2 = (x2(t) - x1(t))^2 + (y2(t) - y1(t))^2 - L2^2; % Enforces constant lenght of second link

% Lagrangian
L = E_kin - E_pot + lambda1*phi1 + lambda2*phi2;

% % Calculate gradients
dLdr = [diff(L,x1(t)); diff(L,y1(t)); diff(L,x2(t)); diff(L,y2(t))]; 
dLdr_dot = [diff(L,x1_dot); diff(L,y1_dot); diff(L,x2_dot); diff(L,y2_dot)];
ddLdr_dotdt = diff(dLdr_dot,t);

% Euler-Lagrange equation
EuLa_eq = ddLdr_dotdt - dLdr; 

% Substitute second derivatives with symbolic variables so we can solve for them
EuLa_eq_sub = subs(EuLa_eq, ...
    [diff(x1(t),t,2), diff(y1(t),t,2), diff(x2(t),t,2), diff(y2(t),t,2)], ...
    [x1dd, y1dd, x2dd, y2dd]);

% Solve for the second derivatives
sol = solve(EuLa_eq_sub == 0, [x1dd, y1dd, x2dd, y2dd]);
x1dd_expr = sol.x1dd;
y1dd_expr = sol.y1dd;
x2dd_expr = sol.x2dd;
y2dd_expr = sol.y2dd;

% Group into vector
r_dot2 = [x1dd_expr; y1dd_expr; x2dd_expr; y2dd_expr];
r1_dot2 = [x1dd_expr; y1dd_expr];
r2_dot2 = [x2dd_expr; y2dd_expr];

% Use the second derivative of the first pendulum link constraint
Cnstr1_eq = r1.'*r1_dot2  +  r1_dot.'*r1_dot;
Cnstr2_eq = s.'*(r2_dot2 - r1_dot2)  +  (r2_dot-r1_dot).'*(r2_dot-r1_dot);

% Solve for the Lagrange multipliers
lambda_sol = solve([Cnstr1_eq == 0, Cnstr2_eq == 0], [lambda1,lambda2]);
lambda1_sol = lambda_sol.lambda1;
lambda2_sol = lambda_sol.lambda2;


% Substitute lambda back into the second order ODE
d2adt = simplify(subs(r_dot2, [lambda1, lambda2], [lambda1_sol, lambda2_sol]));

% Final set of first order ODE's 
dxdt = [r_dot; d2adt];


%%% Create now a .m function to use in a ODE solver after passing anonym
% First substitute the derivatives then the states. This order seems to matter. Matlab...
dxdt_sym = subs(dxdt, ...
   [ diff(x1(t),t), diff(y1(t),t), diff(x2(t),t), diff(y2(t),t)], ...
   [           X1d,           Y1d,           X2d,           Y2d]);

dxdt_sym = subs(dxdt_sym, ...
   [x1(t),  y1(t),  x2(t),  y2(t)], ...
   [   X1,     Y1,     X2,     Y2]);

% Group the variables
x = [X1; Y1; X2; Y2; X1d; Y1d; X2d; Y2d];
theta = [m1; L1; m2; L2; g];
% Write the function
matlabFunction(dxdt_sym, ...
    'File', 'pendulum_double_cartesian_rhs', ...
    'Vars', { x, theta });





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% For a DOUBLE pendulum in ANGLE coordinates with the pivot at the origin
syms ang1(t) ang2(t) m1 m2 L1 L2 g t real
syms a1dd a2dd real % used to inter-store the second derivatives 
syms alpha beta dalpha_dt dbeta_dt real % For the final ODE's writen to .m

% Construct symbolic position, velocity and acceleration vectors
ang1_d = diff(ang1(t), t);
ang2_d = diff(ang2(t), t);

% Define kinetic and potential energy
E_kin = 1/2*m1*ang1_d^2*L1^2  +  1/2*m2*(ang1_d^2*L1^2 + ang2_d^2*L2^2 + 2*L1*L2*ang1_d*ang2_d*cos(ang1(t) - ang2(t)));
E_pot = -m1*g*L1*cos(ang1(t))  +  -m2*g*(L1*cos(ang1(t)) + L2*cos(ang2(t)));

% Lagrangian
L = E_kin - E_pot;

% Compute gradients 
dLdr = gradient(L, [ang1(t); ang2(t)]); % Same as [diff(L,ang1(t)); diff(L,ang2(t))]; 
dLdr_dot = gradient(L, [ang1_d; ang2_d]); % Same as [diff(L,ang1_d); diff(L,ang2_d)];
ddLdr_dotdt = diff(dLdr_dot,t);

% Euler-Lagrange equation
EuLa_eq = ddLdr_dotdt - dLdr; 

% Substitute second derivatives with symbolic variables so we can solve for them
EuLa_eq_sub = subs(EuLa_eq, ...
    [diff(ang1(t),t,2), diff(ang2(t),t,2)], ...
    [a1dd, a2dd]);

% Solve for the second derivatives
sol = solve(EuLa_eq_sub == 0, [a1dd, a2dd]);
ang1_dd = simplify(sol.a1dd);
ang2_dd = simplify(sol.a2dd);

% Final set of first order ODE's 
dxdt = [ang1_d; ang2_d; ang1_dd; ang2_dd];


%%% Create now a .m function to use in a ODE solver after passing anonym
% First substitute the derivatives then the states. This order seems to matter. Matlab...
dxdt_sym = subs(dxdt, ...
   [ diff(ang1(t),t), diff(ang2(t),t)], ...
   [       dalpha_dt,        dbeta_dt]);

dxdt_sym = subs(dxdt_sym, ...
   [ang1(t),  ang2(t)], ...
   [alpha,  beta]);


% Group the variables
x = [alpha; beta; dalpha_dt; dbeta_dt];
theta = [m1; L1; m2; L2; g];

% Write the right-hand side of the ODE as function
matlabFunction(dxdt_sym, ...
    'File', 'pendulum_double_angle_rhs', ...
    'Vars', { x, theta });


%%% Jacobian
J = jacobian(dxdt_sym, x);

% Write the acobain 
matlabFunction(J, ...
    'File', 'pendulum_double_angle_J', ...
    'Vars', {t, x, theta });
