%% beta_max FROM SIMULATION: beta_max vs alpha and beta_max vs N
% Companion to animateLIPM_walkingAngle.m (same LIPM, same local walking-
% frame recursion, same passive mirror law, same initial conditions).
%
% beta_max is FOUND by simulation (coarse scan over beta + bisection) using
% the path-deviation criterion
%     D(N,beta) = |transverse offset of N-th foothold from the straight line|
%     beta_max  = largest beta in [0,90] deg with D <= Dmax
%
% The closed form is overlaid on every plot as a dashed line:
%
%     beta_max = asin( 2*Dmax / ( N*(N+1)*dL0 ) )
%
% where dL0 = per-step step-length growth at beta = 0 (closed form from the
% one-step transition, same as predictedStepGrowth in the decline script):
%     dL0 = -[(I+M) f]_1 = 2 zc tan(alpha) (cosh(omega*tau) - 1)

clear; clc; close all;

%% ---------------- Model parameters ----------------------------------------
p.g          = 9.807;
p.Ts         = 0.005;
p.zModel     = 0.68;
p.stepLength = 0.2;
p.dy_mid     = 0.5;
p.x0         = 0.12;

%% ---------------- Criterion settings --------------------------------------
Dmax      = 1.0;        % corridor half-width [m]
devMetric = 'foothold'; % 'foothold' or 'com'

%% ---------------- Sweep settings ------------------------------------------
alphaRefDeg = 15;
Nref        = 10;

alphaVals   = [1 2 3 5 8 10 12 15 18 20 25 30 35 40];
Nvals       = [2 3 4 5 6 8 10 12 15 20 25 30];

betaScanDeg  = 0:1:90;
bisectTolDeg = 1e-4;

% LaTeX strings for the plots
eqLatex  = '$\beta_{max}=\sin^{-1}\!\left(\frac{2D_{max}}{N(N+1)\,\Delta L_0}\right)$';
dL0Latex = '$\Delta L_0 = 2 z_c \tan\alpha\,(\cosh\omega\tau-1)$';

%% ================= Sweep 1: beta_max vs alpha (N = Nref) ==================
nA = numel(alphaVals);
betaMaxSim_alpha = zeros(1,nA);
betaMaxAna_alpha = zeros(1,nA);
tau_alpha        = zeros(1,nA);
dL0_alpha        = zeros(1,nA);

fprintf('=== beta_max vs alpha   (N = %d, Dmax = %.2f m, metric = %s) ===\n', Nref, Dmax, devMetric);
fprintf('%8s %10s %10s %14s %14s %10s\n', 'alpha', 'tau [s]', 'dL0 [m]', 'beta_max SIM', 'beta_max ANA', 'diff [%]');
for i = 1:nA
    betaMaxSim_alpha(i) = findBetaMaxSim(alphaVals(i), Nref, Dmax, p, devMetric, betaScanDeg, bisectTolDeg);
    [betaMaxAna_alpha(i), tau_alpha(i), dL0_alpha(i)] = betaMaxClosedForm(alphaVals(i), Nref, Dmax, p);
    fprintf('%7.1f%s %10.3f %10.4f %13.4f%s %13.4f%s %10.4f\n', alphaVals(i), char(176), tau_alpha(i), dL0_alpha(i), ...
        betaMaxSim_alpha(i), char(176), betaMaxAna_alpha(i), char(176), ...
        100*(betaMaxSim_alpha(i)-betaMaxAna_alpha(i))/betaMaxAna_alpha(i));
end

%% ================= Sweep 2: beta_max vs N (alpha = alphaRefDeg) ===========
nN = numel(Nvals);
betaMaxSim_N = zeros(1,nN);
betaMaxAna_N = zeros(1,nN);

fprintf('\n=== beta_max vs N   (alpha = %.1f deg, Dmax = %.2f m, metric = %s) ===\n', alphaRefDeg, Dmax, devMetric);
fprintf('%8s %14s %14s %10s\n', 'N', 'beta_max SIM', 'beta_max ANA', 'diff [%]');
for i = 1:nN
    betaMaxSim_N(i) = findBetaMaxSim(alphaRefDeg, Nvals(i), Dmax, p, devMetric, betaScanDeg, bisectTolDeg);
    betaMaxAna_N(i) = betaMaxClosedForm(alphaRefDeg, Nvals(i), Dmax, p);
    fprintf('%8d %13.4f%s %13.4f%s %10.4f\n', Nvals(i), ...
        betaMaxSim_N(i), char(176), betaMaxAna_N(i), char(176), ...
        100*(betaMaxSim_N(i)-betaMaxAna_N(i))/betaMaxAna_N(i));
end
fprintf('\nNote: beta_max = 90 means "unconstrained" (D never exceeds Dmax for any beta in [0,90]).\n');
fprintf('Note: tau differs with alpha because findInitialConditions sizes the stride with gEff = g*cos(alpha).\n');

%% ================= Plot 1: beta_max vs alpha ==============================
figure('Name','beta_max vs alpha (simulated)');
plot(alphaVals, betaMaxAna_alpha, 'w:', 'LineWidth', 2); hold on
plot(alphaVals, betaMaxSim_alpha, 'ro', ...
    'MarkerFaceColor', 'r', 'MarkerSize', 7, 'LineStyle', 'none');

xlabel('\alpha, incline angle (deg)');
ylabel('\beta_{max} (deg)');
title(sprintf('\\beta_{max} vs \\alpha  (N = %d, D_{max} = %.2f m, metric: %s)', Nref, Dmax, devMetric));
legend({['closed form: ' eqLatex], 'simulation'}, 'Interpreter', 'latex', 'Location', 'northeast');
grid on; ylim([0 95]);


%% ================= Plot 2: beta_max vs N ==================================
figure('Name','beta_max vs N (simulated)');
plot(Nvals, betaMaxAna_N, 'w:', 'LineWidth', 2); hold on
plot(Nvals, betaMaxSim_N, 'ro', ...
    'MarkerFaceColor', 'r', 'MarkerSize', 7, 'LineStyle', 'none');

xlabel('N (step horizon)');
ylabel('\beta_{max} (deg)');
title(sprintf('\\beta_{max} vs N  (\\alpha = %.0f%s, D_{max} = %.2f m)', alphaRefDeg, char(176), Dmax));
legend({['closed form: ' eqLatex], 'simulation'}, 'Interpreter', 'latex', 'Location', 'northeast');
grid on; ylim([0 95]);

%% ===================== Local functions =====================================
function bMax = findBetaMaxSim(alphaDeg, N, Dmax, p, metric, betaScanDeg, tolDeg)
    nScan = numel(betaScanDeg);
    D = zeros(1,nScan);
    for k = 1:nScan
        D(k) = simulateDeviation(alphaDeg, betaScanDeg(k), N, p, metric);
    end
    idx = find(D > Dmax, 1, 'first');
    if isempty(idx)
        bMax = betaScanDeg(end);
        return
    end
    if idx == 1
        bMax = betaScanDeg(1);
        return
    end
    lo = betaScanDeg(idx-1);
    hi = betaScanDeg(idx);
    while (hi - lo) > tolDeg
        mid = 0.5*(lo + hi);
        if simulateDeviation(alphaDeg, mid, N, p, metric) > Dmax
            hi = mid;
        else
            lo = mid;
        end
    end
    bMax = 0.5*(lo + hi);
end

function dev = simulateDeviation(alphaDeg, betaDeg, N, p, metric)
    th = deg2rad(alphaDeg);
    bt = deg2rad(betaDeg);

    gEff    = p.g*cos(th);
    gDrive  = p.g*sin(th);
    gDriveS =  gDrive*cos(bt);
    gDriveQ = -gDrive*sin(bt);

    A = [0           1 0           0 0
         gEff/p.zModel 0 0         0 gDriveQ
         0           0 0           1 0
         0           0 gEff/p.zModel 0 gDriveS
         0           0 0           0 0];
    Ad = expm(A*p.Ts);

    [dx0, y0, dy0, tau] = findInitialConditions(p.stepLength, p.dy_mid, p.x0, p.zModel, gEff, p.Ts);
    nSamp = round(tau/p.Ts) + 1;

    state  = [p.x0; dx0; y0; dy0; 1];
    fx     = -p.x0;
    comMax = 0;

    for stp = 1:N
        s = state;
        for k = 1:nSamp-1
            s = Ad*s;
            comMax = max(comMax, abs(s(1) + fx));
        end
        xf  = s(1);  dxf = s(2);  yf = s(3);  dyf = s(4);
        state = [-xf; dxf; -yf; dyf; 1];
        fx = fx + 2*xf;
    end

    switch lower(metric)
        case 'foothold'
            nominalFx = -p.x0*(-1)^N;
            dev = abs(fx - nominalFx);
        case 'com'
            dev = comMax;
        otherwise
            error('devMetric must be ''foothold'' or ''com''.');
    end
end

function [bMaxDeg, tau, dL0] = betaMaxClosedForm(alphaDeg, N, Dmax, p)
    % Closed form:  beta_max = asin( 2*Dmax / ( N*(N+1)*dL0 ) )
    % dL0 = per-step step-length growth at beta = 0, from the one-step map
    %       [y0;dy0] -> M*[y0;dy0] + f  with the mirror law:  dL0 = -[(I+M)f]_1
    th     = deg2rad(alphaDeg);
    gEff   = p.g*cos(th);
    gDrive = p.g*sin(th);
    [~, ~, ~, tau] = findInitialConditions(p.stepLength, p.dy_mid, p.x0, p.zModel, gEff, p.Ts);
    omega  = sqrt(gEff/p.zModel);
    c = cosh(omega*tau);  s = sinh(omega*tau);
    M = [-c, -s/omega; omega*s, c];
    f = [-(gDrive/omega^2)*(c-1); (gDrive/omega)*s];
    v = -((eye(2) + M)*f);
    dL0 = v(1);                                   % = 2 zc tan(alpha) (cosh(omega tau) - 1)
    arg = 2*Dmax / (N*(N+1)*dL0);
    bMaxDeg = asind(min(arg, 1));
end

function [dx0, y0, dy0, tsinglesupport] = findInitialConditions(stepLength, dy_mid, x0, zModel, g, Ts)
    Tc = sqrt(zModel/g);
    y_mid = 0;
    E = -g/(2*zModel)*y_mid^2 + 0.5*dy_mid^2;
    y0 = -stepLength/2;
    dy0 = sqrt(2*(E+g/(2*zModel)*y0^2));
    tsinglesupport = 2*asinh(stepLength/2/(Tc*dy_mid))*Tc;
    tsinglesupport = floor(tsinglesupport/Ts)*Ts;
    tf = tsinglesupport/2;
    dx0 = -x0/Tc * sinh(tf/Tc) / cosh(tf/Tc);
end