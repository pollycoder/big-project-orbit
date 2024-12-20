%-------------------------------------------------------------------%
% Orbit Module for Spacecraft Design Big Project (GPOPS-II)         %
% Endpoint Function                                                 %
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
function phaseout = EBTransferEndpoint(input)
xf = input.phase(1).finalstate;
tf = input.phase(1).finaltime;
phaseout.objective = -xf(7);
end