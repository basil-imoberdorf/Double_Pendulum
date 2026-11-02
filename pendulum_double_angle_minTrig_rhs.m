function dxdt_sym = pendulum_double_angle_minTrig_rhs(in1,in2)
L1 = in2(2,:);
L2 = in2(4,:);
alpha = in1(1,:);
beta = in1(2,:);
dalpha_dt = in1(3,:);
dbeta_dt = in1(4,:);
g = in2(5,:);
m1 = in2(1,:);
m2 = in2(3,:);


delta = alpha - beta;
% Make the minimum necessary trig calls
sa = sin(alpha);
sb = sin(beta);
sd = sin(delta);
cd = cos(delta);

% Now we can rewrite two trig calls with stuff we already computed using some trig identiy
sd2 = 2*sd.*cd;          % sin(2*(alpha-beta))
cd2 = 1 - 2*sd.^2;       % cos(2*(alpha-beta))

da2 = dalpha_dt.^2;
db2 = dbeta_dt.^2;

% Acceleration of the first angle
den1 = L1.*(m1 + m2 - m2.*cd.^2);

a1 = -( g.*(m1+m2).*sa ...
      - g.*m2.*sb.*cd ...
      + L2.*m2.*db2.*sd ...
      + L1.*m2.*da2.*cd.*sd ) ./ den1;


% Acceleration of the second angle
den2 = L2.*( m1 + 0.5*m2 - 0.5*m2.*cd2 );

a2 = ( -0.5*g.*(m1+m2).*sb ...
       + 0.5*g.*(m1+m2).*sin(alpha+delta) ...
       + L1.*(m1+m2).*da2.*sd ...
       + 0.5*L2.*m2.*db2.*sd2 ) ./ den2;


dxdt_sym = [dalpha_dt;dbeta_dt;a1;a2];
end