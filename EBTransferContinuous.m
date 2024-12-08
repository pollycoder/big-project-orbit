%-------------------------------------------------------------------%
% Orbit Module for Spacecraft Design Big Project (GPOPS-II)         %
% Continuous Function - Dynamics                                    %
%-------------------------------------------------------------------%
% Auxdata (International unit):                                     %
% Isp   =   3100            (s)                                     %
% g0    =   9.8             (m/s^2)                                 %
% AU    =   1.49597871e11   (m)                                     %
% muSun =   1.32712440018   (m^3/s^2)                               %
% Tmax  =   0.2             (N)                                     %
%-------------------------------------------------------------------%
% Input (m, s, kg):                                                 %
% r: N x 3                                                          %
% v: N x 3                                                          %
% m: N x 1                                                          %
%-------------------------------------------------------------------%
% Control (N):                                                      %
% u: N x 3                                                          %
%-------------------------------------------------------------------%
% Output (m, s, kg):                                                %
%-------------------------------------------------------------------%
function phaseout = EBTransferContinuous(input)
%-------------------------------------------------------------------%
%----------------------- Unit Transformation -----------------------%
%-------------------------------------------------------------------%
% Unit: au, year, kg
% km2m = 1000;
% year2day = 365.2422;
% day2sec = 86400;
% au2m = input.auxdata.au_;

Isp = input.auxdata.Isp_;
g0 = input.auxdata.g0_;
Tmax = input.auxdata.Tmax_;
muSun = input.auxdata.muSun_;

%-------------------------------------------------------------------%
%------------------- Describe the Dynamic Equation -----------------%
%-------------------------------------------------------------------%
for iphase = 1:length(input.phase)
% State and control
x = input.phase(iphase).state;
u = input.phase(iphase).control;

r = x(:, 1:3);
v = x(:, 4:6);
m = x(:, 7);

num_points = size(r, 1);
u_norm = zeros(num_points, 1);
drdt = v;
dvdt = zeros(num_points, 3);
dmdt = zeros(num_points, 1);

for i = 1:num_points
    r_current = r(i, :);
    m_current = m(i);
    u_current = u(i, :);
    rnorm_current = norm(r_current);
    dvdt(i, :) = -(muSun/rnorm_current^3) .* r_current + (Tmax/m_current) .* u_current;
    dmdt(i) = -(Tmax/Isp/g0) .* norm(u_current);
    u_norm(i) = norm(u_current);
end

phaseout(iphase).dynamics = [drdt, dvdt, dmdt];
phaseout(iphase).path = u_norm;
end

end