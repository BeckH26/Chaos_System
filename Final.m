clear; clc; close all;

% tweak these to change how the sim behaves
L0 = 1.5;        % natural (rest) length of each spring [m]
R = 2;           % ring radius [m]
c = 0.05;        % damping coefficient
nSprings = 3;    % how many anchors around the ring
nMasses = 2;     % how many masses to simulate

% Euler method time settings
h = 0.001;       % time step [s]
T = 100;         % total simulation time [s]
N = round(T/h);  % number of steps

% video stuff
MAKE_VIDEO = true;
VIDEO_FILE = 'spring_sim.mp4';
VIDEO_FPS = 60;
ANIM_STEP = 80;  % speed

% initial conditions x0 y0 vx0 vy0
ICs = [ 1.1,  1.1,  0,    0.0;
        1.0,  1.0,  0,    0.0;
       -1.0,  0.5,  0.0, -0.5];

% one row per mass
K = [1, 1, 1.5;    % stiffnesses mass 1
     1, 1, 1];     % stiffnesses mass 2

% mass
M = [1, 1];

colors = lines(nMasses);

% place anchors evenly around the ring starting at the top
angles = 90 + (0:nSprings-1) * (360/nSprings);
anchors = R * [cosd(angles)', sind(angles)'];

% pack ICs into row 1 of the state matrix
s = zeros(N, 4*nMasses);
for p = 1:nMasses
    idx = (p-1)*4 + 1;
    s(1, idx)   = ICs(p,1);
    s(1, idx+1) = ICs(p,2);
    s(1, idx+2) = ICs(p,3);
    s(1, idx+3) = ICs(p,4);
end

t = (0:N-1)' * h; % time vector

% Euler's method
for n = 1:N-1
    ds = eom(t(n), s(n,:)', M, K, L0, c, anchors, nMasses);
    s(n+1,:) = s(n,:) + h * ds';
end

% unpack solver output and compute energy for each mass
Pos = cell(nMasses, 1);
for p = 1:nMasses
    cols = (p-1)*4 + (1:4);
    Pos{p}.x  = s(:, cols(1));
    Pos{p}.y  = s(:, cols(2));
    Pos{p}.vx = s(:, cols(3));
    Pos{p}.vy = s(:, cols(4));

    PE = zeros(N, 1);
    for i = 1:nSprings
        d = sqrt((Pos{p}.x - anchors(i,1)).^2 + (Pos{p}.y - anchors(i,2)).^2);
        PE = PE + 0.5 * K(p,i) * (d - L0).^2;
    end
    Pos{p}.KE = 0.5 * M(p) * (Pos{p}.vx.^2 + Pos{p}.vy.^2);
    Pos{p}.PE = PE;
    Pos{p}.Etot = Pos{p}.KE + PE;
end

% build legend
legLabels = cell(nMasses, 1);
for p = 1:nMasses
    kRow = sprintf('%.2f ', K(p,:));
    legLabels{p} = sprintf('Mass %d  (x_0=%.1f, y_0=%.1f, m=%.2f, K=[%s])', ...
                            p, ICs(p,1), ICs(p,2), M(p), strtrim(kRow));
end

mStr = ['M=[' sprintf('%.2f ', M) ']'];

% trajectory plot
figure('Name', 'Trajectory', 'Color', 'w', 'Position', [100 100 600 600]);
hold on; axis equal; axis off;

th = linspace(0, 2*pi, 300);
fill(R*cos(th), R*sin(th), [0.95 0.95 0.98], 'EdgeColor', [0.5 0.5 0.6], 'LineWidth', 1.5);

for i = 1:nSprings
    plot(anchors(i,1), anchors(i,2), 'ko', 'MarkerFaceColor', [0.2 0.2 0.8], 'MarkerSize', 10);
    kVals = arrayfun(@(p) K(p,i), 1:nMasses);
    kStr = ['k' num2str(i) '=[' sprintf('%.2f ', kVals) ']'];
    text(anchors(i,1)*1.12, anchors(i,2)*1.12, kStr, ...
         'FontSize', 7, 'HorizontalAlignment', 'center', 'Color', [0.2 0.2 0.8]);
end
anchorH = plot(nan, nan, 'ko', 'MarkerFaceColor', [0.2 0.2 0.8], 'MarkerSize', 10);

trailH = gobjects(nMasses, 1);
for p = 1:nMasses
    x = Pos{p}.x; y = Pos{p}.y;
    n = length(x);
    seg = max(1, floor(n/1000));
    for j = 1:seg:n-1
        plot(x(j:min(j+seg,n)), y(j:min(j+seg,n)), '-', 'Color', colors(p,:), 'LineWidth', 0.6);
    end
    plot(x(1),   y(1),   'g^', 'MarkerFaceColor', 'g', 'MarkerSize', 9);
    plot(x(end), y(end), 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 9);
    trailH(p) = plot(nan, nan, '-', 'Color', colors(p,:), 'LineWidth', 2);
end
startH = plot(nan, nan, 'g^', 'MarkerFaceColor', 'g', 'MarkerSize', 9);
endH = plot(nan, nan, 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 9);

for i = 1:nSprings
    plot([Pos{1}.x(1) anchors(i,1)], [Pos{1}.y(1) anchors(i,2)], '--', ...
         'Color', [0.6 0.6 0.6], 'LineWidth', 1);
end

title(sprintf('Trajectory  Euler h=%.4f  %s  L_0=%.2f  c=%.3f', h, mStr, L0, c), 'FontSize', 12);
legend([anchorH; trailH; startH; endH], ...
       [{'Spring'}; legLabels; {'Start'}; {'End'}], 'Location', 'southeast');

% phase portraits, figure 2
figure('Name', 'Phase Portraits', 'Color', 'w', 'Position', [720 100 800 400]);
for p = 1:nMasses
    n = length(t);
    seg = max(1, floor(n/1000));
    idx = 1:seg:n;
    c_rep = repmat(colors(p,:), numel(idx), 1);
    subplot(1,2,1); hold on;
    scatter(Pos{p}.x(idx), Pos{p}.vx(idx), 2, c_rep, 'filled');
    subplot(1,2,2); hold on;
    scatter(Pos{p}.y(idx), Pos{p}.vy(idx), 2, c_rep, 'filled');
end
subplot(1,2,1);
xlabel('x'); ylabel('v_x'); title('Phase portrait  x–v_x'); grid on; box on;
dh = gobjects(nMasses, 1);
for p = 1:nMasses
    dh(p) = plot(nan, nan, 'o', 'Color', colors(p,:), 'MarkerFaceColor', colors(p,:), 'MarkerSize', 6);
end
legend(dh, legLabels, 'Location', 'best');
subplot(1,2,2);
xlabel('y'); ylabel('v_y'); title('Phase portrait  y–v_y'); grid on; box on;

% energy over time, figure 3
figure('Name', 'Energy', 'Color', 'w', 'Position', [100 720 800 300]);
hold on;
for p = 1:nMasses
    plot(t, Pos{p}.KE,   'b',  'LineWidth', 1,   'DisplayName', sprintf('Kinetic M%d',   p));
    plot(t, Pos{p}.PE,   'r',  'LineWidth', 1,   'DisplayName', sprintf('Potential M%d', p));
    plot(t, Pos{p}.Etot, 'k',  'LineWidth', 1.5, 'DisplayName', sprintf('Total M%d',     p));
end
xlabel('Time [s]'); ylabel('Energy [J]');
title('Energy over time'); legend; grid on; box on;

% heat map figure 4
figure('Name', 'Potential Landscape', 'Color', 'w', 'Position', [920 720 500 450]);
[Xg, Yg] = meshgrid(linspace(-R,R,300), linspace(-R,R,300));
Vg = zeros(size(Xg));
for i = 1:nSprings
    dg = sqrt((Xg - anchors(i,1)).^2 + (Yg - anchors(i,2)).^2);
    Vg = Vg + 0.5 * K(1,i) * (dg - L0).^2;
end
outside = (Xg.^2 + Yg.^2) > R^2;
Vg(outside) = NaN;
contourf(Xg, Yg, Vg, 40, 'LineStyle', 'none');
colormap(hot); colorbar; hold on;
plot(anchors(:,1), anchors(:,2), 'co', 'MarkerFaceColor', 'c', 'MarkerSize', 8);
for p = 1:nMasses
    n = length(t);
    seg = max(1, floor(n/1000));
    plot(Pos{p}.x(1:seg:end), Pos{p}.y(1:seg:end), '.', 'Color', colors(p,:), 'MarkerSize', 1);
end
title('Potential energy landscape (mass 1 K) + trajectories', 'FontSize', 11);
axis equal; xlabel('x'); ylabel('y');

% animation
figAnim = figure('Name', 'Animation', 'Color', 'w', 'Position', [150 150 650 650]);

if MAKE_VIDEO
    vw = VideoWriter(VIDEO_FILE, 'MPEG-4');
    vw.FrameRate = VIDEO_FPS;
    open(vw);
end

nFrames = floor(N / ANIM_STEP);

% set up the static background: ring + anchors
ax = axes('Parent', figAnim);
hold(ax, 'on'); axis(ax, 'equal'); axis(ax, 'off');
axis(ax, [-R-0.3 R+0.3 -R-0.3 R+0.3]);

th = linspace(0, 2*pi, 300);
fill(R*cos(th), R*sin(th), [0.95 0.95 0.98], ...
     'EdgeColor', [0.5 0.5 0.6], 'LineWidth', 1.5, 'Parent', ax);
plot(ax, anchors(:,1), anchors(:,2), 'o', ...
     'MarkerFaceColor', [0.2 0.2 0.8], 'MarkerEdgeColor', 'k', 'MarkerSize', 10);
for i = 1:nSprings
    kVals = arrayfun(@(p) K(p,i), 1:nMasses);
    kStr = ['k' num2str(i) '=[' sprintf('%.2f ', kVals) ']'];
    text(ax, anchors(i,1)*1.12, anchors(i,2)*1.12, kStr, ...
         'FontSize', 7, 'HorizontalAlignment', 'center', 'Color', [0.2 0.2 0.8]);
end

% create handles for the moving parts
animTrail = gobjects(nMasses, 1);
animMass = gobjects(nMasses, 1);
animSpring = gobjects(nMasses, nSprings);
for p = 1:nMasses
    markerSz = 8 + 5 * M(p);
    animTrail(p) = plot(ax, nan, nan, '-', 'Color', colors(p,:), 'LineWidth', 1.2);
    animMass(p) = plot(ax, nan, nan, 'o', 'MarkerFaceColor', colors(p,:), ...
                       'MarkerEdgeColor', 'k', 'MarkerSize', markerSz);
    for sp = 1:nSprings
        lw = 0.5 + 1.5 * (K(p,sp) / max(K(:)));
        animSpring(p,sp) = plot(ax, nan, nan, '-', 'Color', colors(p,:)*0.6, 'LineWidth', lw);
    end
end

lh = gobjects(nMasses+1, 1);
for p = 1:nMasses
    lh(p) = plot(ax, nan, nan, '-o', 'Color', colors(p,:), 'LineWidth', 2, ...
                 'MarkerSize', 7, 'MarkerFaceColor', colors(p,:));
end
lh(nMasses+1) = plot(ax, nan, nan, 'o', 'MarkerFaceColor', [0.2 0.2 0.8], ...
                     'MarkerEdgeColor', 'k', 'MarkerSize', 9);
legend(ax, lh, [legLabels; {'Spring'}], 'Location', 'southeast');

timeTxt = text(ax, -R+0.05, R-0.15, '', 'FontSize', 11, 'FontWeight', 'bold');

% main animation loop
for f = 1:nFrames
    fi = min(f * ANIM_STEP, N);
    for p = 1:nMasses
        x = Pos{p}.x; y = Pos{p}.y;
        set(animTrail(p), 'XData', x(1:fi), 'YData', y(1:fi));
        set(animMass(p),  'XData', x(fi),    'YData', y(fi));
        for sp = 1:nSprings
            set(animSpring(p,sp), 'XData', [x(fi) anchors(sp,1)], ...
                                  'YData', [y(fi) anchors(sp,2)]);
        end
    end
    set(timeTxt, 'String', sprintf('t = %.1f s', t(fi)));
    drawnow;
    if MAKE_VIDEO
        writeVideo(vw, getframe(figAnim));
    end
end

if MAKE_VIDEO; close(vw); end

% equations of motion
function ds = eom(~, s, M, K, L0, c, anchors, nMasses)
    ds = zeros(4*nMasses, 1);
    for p = 1:nMasses
        idx = (p-1)*4 + 1;
        px = s(idx);   py = s(idx+1);
        vx = s(idx+2); vy = s(idx+3);

        Fx = 0; Fy = 0;
        for i = 1:size(anchors,1)
            ax_ = anchors(i,1);
            ay_ = anchors(i,2);
            dx = ax_ - px;
            dy = ay_ - py;
            d = sqrt(dx^2 + dy^2);
            if d > 1e-10
                F = K(p,i) * (d - L0) / d;
                Fx = Fx + F * dx;
                Fy = Fy + F * dy;
            end
        end

        Fx = Fx - c * vx;
        Fy = Fy - c * vy;

        ds(idx) = vx;
        ds(idx+1) = vy;
        ds(idx+2) = Fx / M(p);
        ds(idx+3) = Fy / M(p);
    end
end
