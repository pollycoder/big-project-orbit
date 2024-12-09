%-------------------------------------------------------------------%
% Orbit Module for Spacecraft Design Big Project (GPOPS-II)         %
% Main Function                                                     %
%-------------------------------------------------------------------%
% Auxdata (International unit):                                     %
% Isp   =   3100            (s)                                     %
% g0    =   9.8             (m/s^2)                                 %
% AU    =   1.49597871e11   (m)                                     %
% muSun =   1.32712440018e20(m^3/s^2)                               %
% Tmax  =   0.2             (N)                                     %
%-------------------------------------------------------------------%
% Unit for solving:                                                 %
%   Length: AU                                                      %
%   Time: year                                                      %
%   Mass: % of m0                                                   %
% Target: 317 Roxane                                                %
% Explorer: Tianwen                                                 %
%   Take-off Weight: 8.2t                                           %
% Hall Thruster: single channel HET-450                             %
%   Power level: 105kW                                              %
%   Maximum thrust: 4.6N                                            %
%   Isp: 5100s                                                      %
% Departure: 2035-12-06 20:00:00                                    %
% Arrival: 2041-12-06 20:00:00                                      %  
%-------------------------------------------------------------------%
clc;clear
%-------------------------------------------------------------------%
%---------------------------- Constant -----------------------------%
%-------------------------------------------------------------------%
m0 = 1e5;
mf = 0.8 * m0;
day2sec = 86400;
muSun = 1.32712440018e20;
au = 1.49597871e11;
Isp = 5100;
g0 = 9.8;
Tmax = 4.6;
epoch_jd = juliandate(datetime("2024-10-17 12:00:00"));                                           
epoch0_jd = juliandate(datetime("2000-01-01 12:00:00"));

% Orbit Element 
% Epoch - 2021-07-01
% Epoch0 - 2000-01-01
aAsteroid = 2.285935641764321 * au;
eAsteroid = 0.08612238557068606;
iAsteroid = deg2rad(1.7661971878337);
OmegaAsteroid = deg2rad(151.3313474259361);
omegaAsteroid = deg2rad(186.7199711143541);
MAsteroid_epoch =deg2rad(174.6697087342658);
fAsteroid_epoch = E2f(M2E(MAsteroid_epoch, eAsteroid), eAsteroid);

aEarth = au;
eEarth = 0.0167086;
iEarth = deg2rad(1.578690);
OmegaEarth = deg2rad(174.9);
omegaEarth = deg2rad(288.1);
MEarth_epoch0 = deg2rad(135.27);
fEarth_epoch0 = E2f(M2E(MEarth_epoch0, eEarth), eEarth);




%-------------------------------------------------------------------%
%----------------- Initial and Final Conditions --------------------%
%-------------------------------------------------------------------%
% Initial and final time
t0_datetime = datetime("2035-12-10 20:00:00");
t0_jd = juliandate(t0_datetime);
t0 = t0_jd * day2sec;

tf_datetime = datetime("2049-12-10 20:00:00");
tf_jd = juliandate(tf_datetime);
tf = tf_jd * day2sec;

% Initial state and final state
dt_epoch_f = (tf_jd - epoch0_jd) * day2sec;
fAsteroid_tf = mod(f0dt2ft(fAsteroid_epoch, dt_epoch_f, aAsteroid, eAsteroid, muSun), 2*pi);
coe_tf = [aAsteroid; eAsteroid; iAsteroid; OmegaAsteroid; omegaAsteroid; fAsteroid_tf];
[Rf, Vf] = coe2rv(coe_tf, muSun);

dt_epoch0_0 = (t0_jd - epoch_jd) * day2sec;
fEarth_t0 = mod(f0dt2ft(fEarth_epoch0, dt_epoch0_0, aEarth, eEarth, muSun), 2*pi);
coe_t0 = [aEarth; eEarth; iEarth; OmegaEarth; omegaEarth; fEarth_t0];
[R0, V0] = coe2rv(coe_t0, muSun);

x0 = zeros(1, 7);
x0(1:3) = R0; x0(4:6) = V0; x0(7) = m0;
xf = zeros(1, 7);
xf(1:3) = Rf; xf(4:6) = Vf; xf(7) = mf;                                   % No less than 50% mass remain


%-------------------------------------------------------------------%
%------------------------- Unit Transfer ---------------------------%
%-------------------------------------------------------------------%
tUnit = 365.23 * day2sec;
lUnit = au;
mUnit = m0;
vUnit = lUnit / tUnit;
aUnit = lUnit / tUnit^2;
muUnit = lUnit^3 / tUnit^2;
FUnit = mUnit * aUnit;
m0 = m0 / mUnit;
mf = mf / mUnit;
muSun = muSun / muUnit;
Isp = Isp / tUnit;
g0 = g0 / aUnit;
Tmax = Tmax / FUnit;

% Bounds
rmax = 5;
rmin = -rmax;
vmax = 1e5 / vUnit;
vmin = -vmax;

% Initial and final state/time
x0(1:3) = x0(1:3) ./ lUnit;
x0(4:6) = x0(4:6) ./ vUnit;
x0(7) = x0(7) ./ mUnit;
xf(1:3) = xf(1:3) ./ lUnit;
xf(4:6) = xf(4:6) ./ vUnit;
xf(7) = xf(7) ./ mUnit;
t0 = t0 / tUnit;
tf = tf / tUnit;
coe_t0(1) = coe_t0(1) / lUnit;
coe_tf(1) = coe_tf(1) / lUnit;

%-------------------------------------------------------------------%
%---------------------------- Auxdata ------------------------------%
%-------------------------------------------------------------------%
auxdata.muSun_ = muSun;
auxdata.Isp_ = Isp;
auxdata.g0_ = g0;
auxdata.Tmax_ = Tmax;


%-------------------------------------------------------------------%
%------------------------- Guess: Multi loops ----------------------%
%-------------------------------------------------------------------%
iphase = 1;
NGuess(iphase) = 4;                                                % Number of loops
N = NGuess(iphase);
Node = 1000;                                                                 % Nodes of one loop
guessNode = Node * (N+1);                                                   % Total guess nodes

coe_tf(6) = mod(coe_tf(6), 2*pi) + N*2*pi;

for i=1:guessNode
    coetemp = coe_t0 - (coe_t0-coe_tf) / (guessNode-1) * (i-1);
    coetemp(1) = abs(coetemp(1));
    coetemp(2) = abs(coetemp(2));
    [rtemp, vtemp] = coe2rv(coetemp, muSun);
    mtemp = m0 - (m0-mf) / (guessNode-1) * (i-1);
    guess.phase(iphase).state(i,:) = [rtemp', vtemp', mtemp];
    guess.phase(iphase).control(i,:) = [0.5, 0.5, 0.5];
end
guess.phase(iphase).time = (0:(tf-t0)/(guessNode-1):tf-t0)';


%-------------------------------------------------------------------%
%------------------------------- Bounds ----------------------------%
%-------------------------------------------------------------------%
% Time bounds
bounds.phase(iphase).initialtime.lower = 0;
bounds.phase(iphase).initialtime.upper = 0;
bounds.phase(iphase).finaltime.lower = 0;
bounds.phase(iphase).finaltime.upper = tf-t0;

% State bounds
bounds.phase(iphase).initialstate.lower = x0;
bounds.phase(iphase).initialstate.upper = x0;
bounds.phase(iphase).finalstate.lower = [xf(1:6), xf(7)];
bounds.phase(iphase).finalstate.upper = [xf(1:6), x0(7)];
bounds.phase(iphase).state.lower=[rmin*ones(1,3), vmin*ones(1,3), mf];
bounds.phase(iphase).state.upper=[rmax*ones(1,3), vmax*ones(1,3), m0];

% Control bounds
bounds.phase(iphase).control.lower = -ones(1,3);
bounds.phase(iphase).control.upper = ones(1,3);

% Path constraint
bounds.phase(iphase).path.lower = 0;
bounds.phase(iphase).path.upper = 1;


%-------------------------------------------------------------------%
%-------------------------------- Mesh -----------------------------%
%-------------------------------------------------------------------%
mesh.method = 'hp1';
mesh.tolerance = 1e-11; 
mesh.maxiteration = 10000;
mesh.colpointmin = 4;
mesh.colpointmax = 50;
mesh.phase.colpoints = 10*ones(1,10);
mesh.phase.fraction = 100*ones(1,10);


%-------------------------------------------------------------------%
%--------------------------- Problem Setup -------------------------%
%-------------------------------------------------------------------%
setup.name='EBTransfer-Problem'; 
setup.functions.continuous=@EBTransferContinuous;
setup.functions.endpoint=@EBTransferEndpoint;
setup.auxdata=auxdata;
setup.bounds=bounds;
setup.guess=guess;


%-------------------------------------------------------------------%
%------------------- Solve Problem Using GPOPS2 --------------------%
%-------------------------------------------------------------------%
output = gpops2(setup);
solution = output.result.solution;
result = output.result.objective;
fprintf("J=%f",result);

%%
t = solution.phase.time;
r = solution.phase.state(:, 1:3);
figure
plot3(r(:, 1), r(:, 2), r(:, 3), 'LineWidth', 1.5); hold on
plot3(0, 0, 0, 'k*', 'LineWidth', 3);hold on
plot3(r(1, 1), r(1, 2), r(1, 3), 'g*', 'LineWidth', 3);hold on
text(r(1, 1), r(1, 2), r(1, 3), 'Departure - Earth');hold on
plot3(r(end, 1), r(end, 2), r(end, 3), 'r*', 'LineWidth', 3);hold on
text(r(end, 1), r(end, 2), r(end, 3), 'Arrival - Roxane');
xlabel('x(AU)')
ylabel('y(AU)')
zlabel('z(AU)')
title('Roxane-Earth Low Thrust Trajectory');
axis equal

u = solution.phase.control;
unorm = sqrt(u(:, 1).^ 2 + u(:, 2).^2 + u(:, 3).^2);
figure
plot(t, unorm, 'LineWidth', 1.5);
xlabel('t(Year)')
ylabel('u')
title('Amplitude of Thrust - Normalized')





