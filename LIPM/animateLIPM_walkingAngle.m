%% Simulate LIPM Walking at an Angle Within an Inclined Plane
% Extension of animateLIPM_decline_growth.mlx (which is itself an
% extension of animateLIPM.mlx). Adds a second angle, walkingAngleDeg
% (beta), measured between the intended walking direction and the
% steepest-downhill direction, both measured within the inclined plane.
%
%   declineAngleDeg = 0  -> reduces exactly to animateLIPM.mlx
%   walkingAngleDeg = 0  -> reduces exactly to animateLIPM_decline_growth.mlx
%
% MODELING SUMMARY (full derivation in the accompanying writeup):
% The state-space still has 2 planar coordinates + their velocities, but
% they are now LOCAL, walking-frame coordinates, relative to the stance
% foot, exactly as before:
%   y-row (3/4) = s = position along the WALKING direction
%   x-row (1/2) = q = position TRANSVERSE to the walking direction
% The isotropic restoring term gEff/zModel is a rotation-invariant
% quantity (same coefficient multiplies both s and q), so it is UNCHANGED
% by beta -- only the gravity DRIVE term needs decomposing once beta~=0:
%   gDriveS =  gDrive*cos(beta)   (along walking direction, y/s-row)
%   gDriveQ = -gDrive*sin(beta)   (transverse to walking, x/q-row, NEW)
% A rotation by beta (helper local2world, defined below) converts these
% LOCAL walking-frame coordinates into WORLD (X=lateral, Y=downhill)
% coordinates ONLY where true position on the slope is needed: the
% Z-height formula and the 3-D plot. The recursive LIPM/mirror-law
% arithmetic itself stays in the local frame, exactly as in the original.

animateOn = true;
xTorso = 0.12; % x distance from body to foot
g = 9.807;
Ts = 0.005; %sample time (s)

%% --- Incline + walking-direction parameters ---
declineAngleDeg = 15;   % alpha: magnitude of downhill slope (deg)
walkingAngleDeg = 30;   % beta:  angle of walking direction from steepest-downhill (deg)
                        % beta=0 -> straight downhill; beta=90 -> along the contour
kx = 0;                        % plane has no tilt in the world-X direction
ky = -tand(declineAngleDeg);   % world tilt is always along true downhill (Y), regardless of beta

theta = deg2rad(declineAngleDeg);
beta  = deg2rad(walkingAngleDeg);

gEff   = g*cos(theta);   % normal-to-plane gravity -> sets pendulum stiffness gEff/z (unchanged by beta)
gDrive = g*sin(theta);   % total in-plane (downhill) gravity, magnitude fixed by alpha alone

% Decompose the fixed-magnitude gDrive into the walking-frame axes.
% sqrt(gDriveS^2+gDriveQ^2) == gDrive for ANY beta -- beta only decides
% how that fixed forcing splits between "drives me forward" and "pushes
% me sideways" (see writeup section 3 for the geometric derivation).
gDriveS =  gDrive*cos(beta);   % along walking direction (added to the y/s-row)
gDriveQ = -gDrive*sin(beta);   % transverse to walking direction (added to the x/q-row)

%% Robot / model heights
zRobot = 0.78;
zModel = 0.68;

%% Swing
swingHeight = 0.1;

%% Initial conditions
% findInitialConditions still sizes a nominal stride using gEff only
% (unaffected by beta), exactly as in the original decline model.
useParametricICs = true;

if useParametricICs
    stepLength = 0.2;
    dy_mid = 0.5;
    x0 = 0.12;
    [dx0, y0, dy0, tSingleSupport] = findInitialConditions(stepLength, dy_mid, x0, zModel, gEff, Ts);
else
    tSingleSupport = 0.37;
    x0 = 0.12;
    dx0 = -0.4501;
    y0 = -0.1;
    dy0 = 0.3845;
end

%% State-space model (LOCAL walking-frame coordinates)
% Row 1/2: q (transverse-to-walking, local "x") <- NOW carries gDriveQ
% Row 3/4: s (along-walking, local "y")          <- carries gDriveS
A = [0           1 0           0 0
    gEff/zModel  0 0           0 gDriveQ
    0            0 0           1 0
    0            0 gEff/zModel 0 gDriveS
    0            0 0           0 0];
B = [0 0
     1 0
     0 0
     0 1
     0 0];
C = [1 0 0 0 0
     0 0 1 0 0];
D = [0 0; 0 0];

lipm  = ss(A,B,C,D);
lipmD = c2d(lipm, Ts);

Ad = lipmD.A;
Bd = lipmD.B;
Cd = lipmD.C;
Dd = lipmD.D;

%% Initial state
state0 = [x0; dx0; y0; dy0; 1];

%% Closed-form per-step growth prediction (both axes)
omega  = sqrt(gEff/zModel);
tau    = tSingleSupport;
coshWT = cosh(omega*tau);
sinhWT = sinh(omega*tau);
M = [-coshWT, -sinhWT/omega; omega*sinhWT, coshWT];

fS = [-(gDriveS/omega^2)*(coshWT-1); (gDriveS/omega)*sinhWT];
fQ = [-(gDriveQ/omega^2)*(coshWT-1); (gDriveQ/omega)*sinhWT];

predictedGrowthS = -((eye(2)+M)*fS); predictedGrowthS = predictedGrowthS(1); % forward growth/step
predictedGrowthQ = -((eye(2)+M)*fQ); predictedGrowthQ = predictedGrowthQ(1); % transverse drift/step

fprintf('Predicted forward    step-length growth: %.4f m/step\n', predictedGrowthS);
fprintf('Predicted transverse drift growth:       %.4f m/step\n', predictedGrowthQ);
fprintf('Total in-plane growth magnitude:         %.4f m/step (invariant to beta)\n', ...
    hypot(predictedGrowthS,predictedGrowthQ));

%% Stability check: beta_max via a leg-reach (kinematic) limit
% *** MINIMUM ADDITIONAL ASSUMPTION -- the bare LIPM has none of this ***
% The point-foot LIPM + mirror law has NO intrinsic beta-dependent
% instability: the required friction ratio mu = tan(alpha) is IDENTICAL
% for every beta (because the total drive magnitude is fixed by alpha
% alone -- verified numerically in the writeup). A finite beta_max
% therefore requires an extra, explicitly-flagged kinematic assumption:
% the swing leg cannot be placed farther than rMax from the hip. zRobot
% is used as a stand-in for leg length since the model has no dedicated
% leg-length parameter.
legLength = zRobot;                        % ASSUMPTION: leg length ~= standing height
rMax = sqrt(legLength^2 - zModel^2);        % max single-side horizontal foot reach
dL0  = 2*zModel*tan(theta)*(coshWT-1);      % total per-step growth magnitude (beta-independent)

if rMax < dL0
    betaMaxDeg = asind(rMax/dL0);
else
    betaMaxDeg = 90; % leg reach exceeds even the beta=90 transverse growth: unconstrained by this criterion
end
fprintf('\nLeg-reach transverse limit rMax = %.4f m\n', rMax);
fprintf('Maximum feasible walking angle (this criterion) = %.2f deg\n', betaMaxDeg);
if walkingAngleDeg <= betaMaxDeg
    fprintf('Current walkingAngleDeg = %.1f deg -> FEASIBLE under this criterion\n\n', walkingAngleDeg);
else
    fprintf('Current walkingAngleDeg = %.1f deg -> EXCEEDS this criterion\n\n', walkingAngleDeg);
end

%% u0
u0 = [0; 0];

%% modelVars
modelVars.g = g;
modelVars.Ts = Ts;
modelVars.zModel = zModel;
modelVars.kx = kx;
modelVars.ky = ky;
modelVars.beta = beta;             % NEW: needed by stanceSim's world-frame conversion
modelVars.Ad = Ad;
modelVars.Bd = Bd;
modelVars.Cd = Cd;
modelVars.Dd = Dd;
modelVars.tSingleSupport = tSingleSupport;

%% simStates
simStates = struct('bodyPos',[],'footPos',[],'timeVec',[]);

%% Initial feet -- placed PERPENDICULAR TO THE WALKING DIRECTION (rotated
% by beta from the start), not fixed to world-X. This is what makes the
% robot actually walk diagonally rather than just having its plane/camera
% rotated.
[Xw,Yw] = local2world(-xTorso, 0, beta);
simStates.footPos(:,end+1) = [Xw; Yw; 0]; % left foot
[Xw,Yw] = local2world( xTorso, 0, beta);
simStates.footPos(:,end+1) = [Xw; Yw; 0]; % right foot

numSteps = 10; % enough steps to see the (linear) growth clearly

%% Figure setup
robotpos0 = [0; 0; zRobot];

hFig = figure;
hFig.Visible = 'on'; % for livescript
hFig.Units = 'Normalized';
hFig.OuterPosition = [0 0 1 1];
hAx  = axes(hFig);

hold(hAx, 'on');
hLeg  = plot3(hAx,[0,robotpos0(1)],[0,robotpos0(2)],[0,robotpos0(3)],...
                '-o','MarkerIndices',2,'LineWidth',2,'MarkerSize',8); % com + leg
hCom = plot3(hAx,robotpos0(1),robotpos0(2),robotpos0(3)); % com trajectory
hLineFootHold = plot3(hAx, simStates.footPos(1,:), simStates.footPos(2,:), simStates.footPos(3,:)); % foothold traj
hLineLeftFoot = plot3(hAx, simStates.footPos(1,1), simStates.footPos(2,1), simStates.footPos(3,1));
hLineRightFoot = plot3(hAx, simStates.footPos(1,2), simStates.footPos(2,2), simStates.footPos(3,2));
hold(hAx, 'off');

grid(hAx, 'on');
view(hAx,3)
axis equal
% Axis limits now must accommodate growth in BOTH world X and Y (since
% walking direction can point anywhere in-plane), sized off dL0, the
% beta-independent total per-step growth magnitude computed above.
halfSpan = 3*zRobot + numSteps*dL0*1.3;
hAx.XLim = [-halfSpan halfSpan];
hAx.YLim = [-halfSpan halfSpan];

% Draw the inclined ground plane for visual reference (always tilted
% along true downhill/world-Y, independent of beta)
hold(hAx, 'on');
[Xg, Yg] = meshgrid(linspace(hAx.XLim(1), hAx.XLim(2), 2), linspace(hAx.YLim(1), hAx.YLim(2), 2));
Zg = kx*Xg + ky*Yg;
surf(hAx, Xg, Yg, Zg, 'FaceColor', [0.8 0.8 0.8], 'FaceAlpha', 0.25, ...
     'EdgeColor', 'none', 'HandleVisibility', 'off');

% Reference arrows: steepest-downhill direction (blue, dashed) and the
% intended walking direction (red), both from the origin.
arrowLen = 1.5*zRobot;
plot3(hAx, [0 0], [0 arrowLen], [0 0], 'b--', 'LineWidth', 1.5, 'HandleVisibility','off');
[Xa,Ya] = local2world(0, arrowLen, beta);
plot3(hAx, [0 Xa], [0 Ya], [0 0], 'r-', 'LineWidth', 2, 'HandleVisibility','off');
text(hAx, 0, arrowLen, 0.05, 'downhill', 'Color','b');
text(hAx, Xa, Ya, 0.05, sprintf('walking dir (\\beta=%.0f%s)', walkingAngleDeg, char(176)), 'Color','r');
hold(hAx, 'off');

legend(hAx, 'Virtual Leg and COM of LIP','Trajectory of COM', 'Foot Hold', 'Left Foot', 'Right Foot', ...
       'Location', 'northeast');

%% stepinfos / stepRecord
stepinfos = cell(1,numSteps+2);
% columns: [step, worldFootholdX, worldFootholdY, forwardStepLen, transverseStepLen]
stepRecord = zeros(numSteps,5);

%% Move COM from standing to over one foot (local transverse shift, rotated to world)
[Xr,Yr] = local2world(0.67*x0, 0, beta);
walkRdyPoints = [0      Xr;
                 0      Yr;
                 zRobot zModel];
% Note: COM is moving partly over the (local) left/first foot.

timepoints = [0 1]; % Hard-coded 1 second to move down
timevec = timepoints(1):Ts:timepoints(end)-Ts;
[q, qd] = cubicpolytraj(walkRdyPoints, timepoints, timevec);
simStates.bodyPos = [simStates.bodyPos q];
simStates.timeVec = [simStates.timeVec timevec];

% save data
fhold_x = simStates.footPos(1,1);
fhold_y = simStates.footPos(2,1);
stepinfo.index = 1;
bodystate = [q(1,:); qd(1,:); q(2,:); qd(2,:); q(3,:); qd(3,:)]; % pos wrt world
relstate = bodystate - [fhold_x; 0; fhold_y; 0; 0; 0]; % save pos wrt foot
stepinfo.state = relstate;
stepinfo.timevec = timevec;
stepinfo.mode = 'doublesupport';
stepinfo.footplant = [fhold_x; fhold_y; 0];
stepinfos{1} = stepinfo;

% update graphics
for idx = 1:size(simStates.bodyPos,2)
    points = [zeros(3,1) simStates.bodyPos(:,idx)];
    updateLine(hLeg, points);
    appendLine(hCom, simStates.bodyPos(:,idx));
    if animateOn
        pause(0.005);
    end
end

%% First half-step: local advance sAdvance along the walking direction
sAdvance = 0.1; % same magnitude as the original script's literal 0.1
[Xr2,Yr2] = local2world(0, sAdvance, beta);
walkRdyPoints(:,end+1) = [Xr2; Yr2; kx*Xr2 + ky*Yr2 + zModel];

% new foothold: LOCAL construction (matches the mirror-law bookkeeping
% used throughout the rest of the script), rotated to world only for storage/plotting
fhold_x_local = -x0;
fhold_y_local = -y0 + sAdvance;
[Xw,Yw] = local2world(fhold_x_local, fhold_y_local, beta);
simStates.footPos(:,end+1) = [Xw; Yw; kx*Xw + ky*Yw]; % left foot (world)

% make trajectory for the stance foot to walk-initial state
rdyPos1 = walkRdyPoints(:,2);
rdypos2 = walkRdyPoints(:,3);
rdyvel1 = [0; 0; 0];
[dXv,dYv] = local2world(dx0, dy0, beta); % rotate the IC velocity too (linear map commutes with d/dt)
rdyvel2 = [dXv; dYv; 0];
waypoints = [rdyPos1 rdypos2];
velpoints = [rdyvel1 rdyvel2];

timepoints = [1 1.5]; % Hard-coded 0.5 seconds to take a half step
timevec = timepoints(1):Ts:timepoints(end)-Ts;
[q, qd, qdd] = cubicpolytraj(waypoints, timepoints, timevec, ...
                  'VelocityBoundaryCondition', velpoints);
simStates.bodyPos = [simStates.bodyPos q];
simStates.timeVec = [simStates.timeVec timevec];

% make trajectory for the swing foot
swingfootpos0 = simStates.footPos(:,1);
swingfootpos1 = simStates.footPos(:,3);
[qswing, dqswing, ddqswing] = getSwingFootTraj(swingfootpos0, swingfootpos1, swingHeight, ...
                                timepoints(1), timepoints(end),Ts);

% save data
stepinfo.index = 2;
bodystate = [q(1,:); qd(1,:); q(2,:); qd(2,:); q(3,:); qd(3,:)];
relstate = bodystate - [simStates.footPos(1,2); 0; simStates.footPos(2,2); 0; 0; 0];
stepinfo.state = relstate;
stepinfo.timevec = timevec;
stepinfo.mode = 'singlesupportright';
stepinfo.swing = qswing;
stepinfo.footplant = [simStates.footPos(1,2); simStates.footPos(2,2); 0];
stepinfos{2} = stepinfo;

% update graphics
for idx = 1:size(q,2)
    points = [zeros(3,1) q(:,idx)];
    updateLine(hLeg, points);
    appendLine(hCom, q(:,idx));
    appendLine(hLineLeftFoot, qswing(:,idx));
    if animateOn
        pause(0.005);
    end
end
appendLine(hLineFootHold, [Xw; Yw; kx*Xw+ky*Yw])

%% Main stepping loop (passive mirror law, LOCAL walking-frame recursion)
fhold_x1 = fhold_x_local;
fhold_y1 = fhold_y_local;

for stp = 1:numSteps
    fhold_x0 = fhold_x1; % LOCAL
    fhold_y0 = fhold_y1; % LOCAL

    [states, timevec, simStates] = stanceSim(fhold_x0, fhold_y0, state0, u0, modelVars, simStates);
    [fholdX0_w, fholdY0_w] = local2world(fhold_x0, fhold_y0, beta); % world, for plotting/logging only

    index = stp + 2;
    stepinfo.index = index;
    [XabsW, YabsW] = local2world(states(1,:) + fhold_x0, states(3,:) + fhold_y0, beta);
    zVec  = kx*XabsW + ky*YabsW + zModel;
    [dXabsW, dYabsW] = local2world(states(2,:), states(4,:), beta);
    dzVec = kx*dXabsW + ky*dYabsW;
    savestate = [states(1:4,:); zVec; dzVec];
    stepinfo.state = savestate;
    stepinfo.footplant = [fholdX0_w; fholdY0_w; kx*fholdX0_w + ky*fholdY0_w];
    stepinfo.timevec = timevec;
    if rem(stp,2) > 0
        stepinfo.mode = 'singlesupportleft';
    else
        stepinfo.mode = 'singlesupportright';
    end

    % change leg according to the passive mirror law (LOCAL frame)
    [state0, fhold_x1, fhold_y1] = changeLeg(states(:,end), fhold_x0, fhold_y0);
    [fholdX1_w, fholdY1_w] = local2world(fhold_x1, fhold_y1, beta);
    simStates.footPos(:,end+1) = [fholdX1_w; fholdY1_w; kx*fholdX1_w + ky*fholdY1_w];

    % report and record: LOCAL step-length components (forward + transverse)
    stepLenS = fhold_y1 - fhold_y0; % along walking direction
    stepLenQ = fhold_x1 - fhold_x0; % transverse drift
    stepRecord(stp,:) = [stp, fholdX1_w, fholdY1_w, stepLenS, stepLenQ];
    fprintf('Step %2d: world foothold = (%6.3f, %6.3f) m   forward step = %5.3f m   transverse step = %6.3f m\n', ...
        stp, fholdX1_w, fholdY1_w, stepLenS, stepLenQ);

    % define swing foot trajectory with new leg position
    [q, dq, ddq] = getSwingFootTraj(simStates.footPos(:,end-2), simStates.footPos(:,end), ...
        swingHeight,timevec(1), timevec(end), Ts);
    stepinfo.swing = q;

    % simulate stance
    for idx = 1:size(states,2)
        state1 = states(:,idx);
        meas1 = Cd*state1 + Dd*u0;

        [Xw, Yw] = local2world(fhold_x0 + meas1(1), fhold_y0 + meas1(2), beta);
        bodypos = [Xw; Yw; kx*Xw + ky*Yw + zModel];

        % animate stance leg and COM traj
        points = [fholdX0_w, bodypos(1);
                  fholdY0_w, bodypos(2);
                  kx*fholdX0_w + ky*fholdY0_w, bodypos(3)];
        updateLine(hLeg, points)
        appendLine(hCom, bodypos);

        % animate swing leg
        if rem(stp,2) > 0
            appendLine(hLineRightFoot, q(:,idx))
        else
            appendLine(hLineLeftFoot, q(:,idx))
        end
        if animateOn
            pause(0.005);
        end
    end
    % update graphics: change leg
    points = [fholdX1_w, bodypos(1);
              fholdY1_w, bodypos(2);
              kx*fholdX1_w + ky*fholdY1_w, bodypos(3)];
    updateLine(hLeg, points)
    appendLine(hLineFootHold, [fholdX1_w; fholdY1_w; kx*fholdX1_w + ky*fholdY1_w])
    if animateOn
        pause(0.005);
    end

    % save the data
    stepinfos{stp+2} = stepinfo;
end

% Summary table of every foothold position (world) and step length (local components)
stepTable = array2table(stepRecord, ...
    'VariableNames', {'Step','FootholdX_world','FootholdY_world','ForwardStep_local','TransverseStep_local'});
disp(stepTable);

% Compare the closed-form prediction against what was actually simulated
observedGrowthS = mean(diff(stepRecord(:,4)));
observedGrowthQ = mean(diff(stepRecord(:,5)));
fprintf('\nPredicted forward   growth/step: %.4f m   Observed: %.4f m\n', predictedGrowthS, observedGrowthS);
fprintf('Predicted transverse growth/step: %.4f m   Observed: %.4f m\n', predictedGrowthQ, observedGrowthQ);

%% ================= Validation: limiting cases =================
fprintf('\n--- Validation ---\n');
fprintf('beta=0 check: gDriveQ = %.6g (should be 0 -> exactly recovers animateLIPM_decline_growth.mlx)\n', gDrive*sin(0));
fprintf('alpha=0 check: gDrive = %.6g (should be 0 for any beta -> exactly recovers animateLIPM.mlx)\n', g*sin(0));
fprintf('beta=90 check: cos(90deg) = %.6g (forward drive vanishes; all growth becomes transverse)\n', cosd(90));

%% ================= Beta sweep (Sections 12) =================
betaSweepDeg = 0:1:90;
nB = numel(betaSweepDeg);
dS_sweep = zeros(1,nB); dQ_sweep = zeros(1,nB);
gS_sweep = zeros(1,nB); gQ_sweep = zeros(1,nB);
for i = 1:nB
    [dS_sweep(i), dQ_sweep(i), gS_sweep(i), gQ_sweep(i)] = ...
        stepGrowthVsAngle(declineAngleDeg, betaSweepDeg(i), zModel, stepLength, dy_mid, x0, g, Ts);
end

figure('Name','Gravity drive components vs walking angle');
plot(betaSweepDeg, gS_sweep, 'r-', 'LineWidth',2); hold on
plot(betaSweepDeg, gQ_sweep, 'b-', 'LineWidth',2);
xlabel('\beta (deg)'); ylabel('drive acceleration (m/s^2)');
legend('gDrive_S (forward)','gDrive_Q (transverse)','Location','best');
title(sprintf('Gravity-drive decomposition, \\alpha=%.0f%s',declineAngleDeg,char(176))); grid on

figure('Name','Step-length growth vs walking angle');
plot(betaSweepDeg, dS_sweep, 'r-', 'LineWidth',2); hold on
plot(betaSweepDeg, dQ_sweep, 'b-', 'LineWidth',2);
plot(betaSweepDeg, hypot(dS_sweep,dQ_sweep), 'k--', 'LineWidth',1.5);
xlabel('\beta (deg)'); ylabel('growth per step (m)');
legend('\DeltaS (forward)','\DeltaQ (transverse)','magnitude (invariant)','Location','best');
title(sprintf('Per-step growth decomposition, \\alpha=%.0f%s',declineAngleDeg,char(176))); grid on

%% ================= beta_max(alpha) sweep (Section 13) =================
alphaSweepDeg = 5:1:30;
betaMaxSweep = zeros(size(alphaSweepDeg));
for i = 1:numel(alphaSweepDeg)
    [~,~,~,~,dL0_i] = stepGrowthVsAngle(alphaSweepDeg(i), 0, zModel, stepLength, dy_mid, x0, g, Ts); %#ok<ASGLU>
    if rMax < dL0_i
        betaMaxSweep(i) = asind(rMax/dL0_i);
    else
        betaMaxSweep(i) = 90;
    end
end

figure('Name','beta_max vs incline angle');
plot(alphaSweepDeg, betaMaxSweep, 'm-', 'LineWidth',2);
xlabel('\alpha, decline angle (deg)'); ylabel('\beta_{max} (deg)');
title(sprintf('Maximum feasible walking angle vs incline (rMax = %.3f m leg-reach limit)', rMax));
grid on; ylim([0 95]);

%% ===================== Local functions =====================
function [Xw, Yw] = local2world(xLocal, yLocal, beta)
    % Rotates LOCAL walking-frame coordinates (xLocal = transverse "q",
    % yLocal = along-walking "s") into WORLD coordinates (X = world
    % lateral, Y = world steepest-downhill). At beta=0 this is the
    % identity map, exactly recovering the original decline model.
    Xw = yLocal.*sin(beta) + xLocal.*cos(beta);
    Yw = yLocal.*cos(beta) - xLocal.*sin(beta);
end

function [dS, dQ, gDriveS, gDriveQ, dL0] = stepGrowthVsAngle(alphaDeg, betaDeg, zModel, stepLength, dy_mid, x0, g, Ts)
    % Closed-form forward/transverse per-step growth at arbitrary alpha, beta.
    th = deg2rad(alphaDeg); bt = deg2rad(betaDeg);
    gEffL = g*cos(th); gDriveL = g*sin(th);
    gDriveS = gDriveL*cos(bt);
    gDriveQ = -gDriveL*sin(bt);
    [~, ~, ~, tau] = findInitialConditions(stepLength, dy_mid, x0, zModel, gEffL, Ts);
    omegaL = sqrt(gEffL/zModel);
    Cw = cosh(omegaL*tau);
    dS = 2*zModel*(gDriveS/gEffL)*(Cw-1);
    dQ = 2*zModel*(gDriveQ/gEffL)*(Cw-1);
    dL0 = 2*zModel*tan(th)*(Cw-1);
end

function [dx0, y0, dy0, tsinglesupport] = findInitialConditions(stepLength, dy_mid, x0, zModel, g, Ts)

    Tc = sqrt(zModel/g);
    % Desired midstance y state
    y_mid = 0;
    % Corresponding orbital energy is
    E = -g/(2*zModel)*y_mid^2 + 0.5*dy_mid^2;

    y0 = -stepLength/2;
    % finding dy0 from midstance energy level
    dy0 = sqrt(2*(E+g/(2*zModel)*y0^2));

    % using relationship between final body state and initial body state,
    % we can find time it will take to reach midstance given final velocity
    % (dy = dy_mid) and final position (which is y = 0 at midstance)
    tsinglesupport = 2*asinh(stepLength/2/(Tc*dy_mid))*Tc;
    tsinglesupport = floor(tsinglesupport/Ts)*Ts; % should be divisible by sample time

    % also using relationship between final body state and initial body
    % state. We already know time it will take to reach midstance.
    tf = tsinglesupport/2;
    dx0 = -x0/Tc * sinh(tf/Tc) / cosh(tf/Tc);
end

function [states, timeVec, simStates] = stanceSim(fhold_x, fhold_y, state0, u0, modelVars, simStates)
    % fhold_x, fhold_y: current stance foot position in LOCAL
    % (walking-frame) coordinates -- the same frame as the LIPM state.

    % Unpack variables
    Ts = modelVars.Ts;
    Ad = modelVars.Ad;
    Bd = modelVars.Bd;
    kx = modelVars.kx;
    ky = modelVars.ky;
    beta = modelVars.beta;

    % Define time vector
    tInitial = simStates.timeVec(end) + modelVars.Ts; % don't let time overlap
    tFinal   = tInitial + modelVars.tSingleSupport;
    timeVec  = tInitial:Ts:tFinal;

    % Simulation loop
    nSteps = numel(timeVec);
    states = zeros(size(state0,1),nSteps);
    states(:,1) = state0;
    for idx = 1:nSteps-1
        states(:,idx+1) = Ad*states(:,idx) + Bd*u0;
    end

    % Update the simStates structure (WORLD frame, via rotation by beta)
    simStates.timeVec = [simStates.timeVec timeVec(1:end-1)];
    xAbsLocal = states(1,2:end) + fhold_x;
    yAbsLocal = states(3,2:end) + fhold_y;
    [XabsW, YabsW] = local2world(xAbsLocal, yAbsLocal, beta);
    zAbs = kx*XabsW + ky*YabsW + modelVars.zModel; % Z must use WORLD Y
    newStates = [XabsW; YabsW; zAbs];
    simStates.bodyPos = [simStates.bodyPos newStates];
end

function [state0, fhold_x, fhold_y] = changeLeg(state1, fhold_x0, fhold_y0)
    % Passive control law (mirroring) -- no braking. Expressed directly
    % in LOCAL walking-frame coordinates. This is an algebraically
    % EQUIVALENT refactor of the original bodyPos-based formula
    % (fhold_new = bodyPos(end) - x0, with bodyPos(end) = fhold_old + xf,
    % x0 = -xf  =>  fhold_new = fhold_old + 2*xf); written this way it no
    % longer needs to assume bodyPos is expressed in the same frame as
    % the state, which matters once a world-frame rotation is involved.
    xf  = state1(1);
    dxf = state1(2);
    yf  = state1(3);
    dyf = state1(4);

    x0 = -xf;
    y0 = -yf;
    dx0 = dxf;
    dy0 = dyf;
    state0 = [x0; dx0; y0; dy0; 1];

    fhold_x = fhold_x0 + 2*xf; % local transverse foothold update
    fhold_y = fhold_y0 + 2*yf; % local forward foothold update
end

function appendLine(gHandle, points)
    gHandle.XData(end+1) = points(1);
    gHandle.YData(end+1) = points(2);
    gHandle.ZData(end+1) = points(3);
end

function updateLine(gHandle, points)
    gHandle.XData = points(1,:);
    gHandle.YData = points(2,:);
    gHandle.ZData = points(3,:);
end

function xi = findFootHold(xd, vd, vi, mvar)
    g = mvar.g; zModel = mvar.zModel; tSingleSupport = mvar.tSingleSupport;
    a = 1; % weight for distance
    b = 1; % weight for velocity

    Tc = sqrt(zModel/g);
    CT = cosh(tSingleSupport/Tc);
    ST = sinh(tSingleSupport/Tc);
    DT = a*CT^2 + b*(ST/Tc)^2;

    xi = (a*CT*(xd - ST*Tc*vd) + b*ST/Tc*(vd-CT*vi))/DT;
end
