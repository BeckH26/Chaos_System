clear; clc; close all;

%  paramaters
m = 1.0; % mass
k = 1.0; % spring stiffness [N/m]
L0 = 1.5; % natural (rest) length of each spring [m]
R = 2.0; % ring radius = distance from centre to each anchor [m]
c = 0.05; % damping coefficient [N·s/m]

x0 = 0.3; % initial x displacement
y0 = 0.8; % initial y displacement
vx0 = 0.5; % initial x velocity
vy0 = -0.2; % initial y velocity

tspan = [0 200];

angles  = [90, 210, 330]; % angles where the springs are
anchors = R * [cosd(angles)', sind(angles)']; % 3x2 matrix

% Equation
ode = @(t, s) eom(t, s, m, k, L0, c, anchors);

%  Integrate
opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-11, 'MaxStep', 0.05);
ic   = [x0; y0; vx0; vy0];
[t, S] = ode45(ode, tspan, ic, opts);

x  = S(:,1);
y  = S(:,2);
vx = S(:,3);
vy = S(:,4);

PE = zeros(size(x));
for i = 1:3
    d   = sqrt((x - anchors(i,1)).^2 + (y - anchors(i,2)).^2);
    PE  = PE + 0.5 * k * (d - L0).^2;
end
KE    = 0.5 * m * (vx.^2 + vy.^2);
Etot  = KE + PE;

% trajecrory in the ring
figure('Name','Trajectory','Color','w','Position',[100 100 600 600]);
hold on; axis equal; axis off;

% Draw ring
th   = linspace(0, 2*pi, 300);
fill(R*cos(th), R*sin(th), [0.95 0.95 0.98], 'EdgeColor',[0.5 0.5 0.6], ...
     'LineWidth', 1.5);

% Draw spring anchors
for i = 1:3
    plot(anchors(i,1), anchors(i,2), 'ko', 'MarkerFaceColor',[0.2 0.2 0.8], ...
         'MarkerSize', 10);
end

% Colour trajectory by time 
n   = length(t);
seg = max(1, floor(n/1000));   % plot ~1000 segments
idx = 1:seg:n;
cmap = jet(length(idx));
for j = 1:length(idx)-1
    i1 = idx(j); i2 = idx(j+1);
    plot(x(i1:i2), y(i1:i2), '-', 'Color', cmap(j,:), 'LineWidth', 0.6);
end

% Mark start
plot(x(1), y(1), 'g^', 'MarkerFaceColor','g', 'MarkerSize', 9);
plot(x(end), y(end), 'rs', 'MarkerFaceColor','r', 'MarkerSize', 9);

% Spring lines at t=0
for i = 1:3
    plot([x(1) anchors(i,1)], [y(1) anchors(i,2)], '--', ...
         'Color',[0.6 0.6 0.6], 'LineWidth', 1);
end

title(sprintf('Trajectory   k=%.2f  L_0=%.2f  c=%.3f', k, L0, c), ...
      'FontSize', 13);
legend({'','Anchor','','Trajectory','Start','End'}, 'Location','southeast');
colormap(jet); cb = colorbar; cb.Label.String = 'Time progression →';
clim([0 1]);

% phase portraits
figure('Name','Phase Portraits','Color','w','Position',[720 100 800 400]);

subplot(1,2,1);
scatter(x(1:seg:end), vx(1:seg:end), 2, t(1:seg:end), 'filled');
colormap(jet); colorbar;
xlabel('x'); ylabel('v_x');
title('Phase portrait  x–v_x'); grid on; box on;

subplot(1,2,2);
scatter(y(1:seg:end), vy(1:seg:end), 2, t(1:seg:end), 'filled');
colormap(jet); colorbar;
xlabel('y'); ylabel('v_y');
title('Phase portrait  y–v_y'); grid on; box on;

%  energy vs time
figure('Name','Energy','Color','w','Position',[100 720 800 300]);
plot(t, KE,   'b', 'LineWidth', 1,   'DisplayName','Kinetic');
hold on;
plot(t, PE,   'r', 'LineWidth', 1,   'DisplayName','Potential');
plot(t, Etot, 'k', 'LineWidth', 1.5, 'DisplayName','Total');
xlabel('Time [s]'); ylabel('Energy [J]');
title('Energy over time'); legend; grid on; box on;


%  heat map
figure('Name','Potential Landscape','Color','w','Position',[920 720 500 450]);
[Xg, Yg] = meshgrid(linspace(-R,R,300), linspace(-R,R,300));
Vg       = zeros(size(Xg));
for i = 1:3
    dg = sqrt((Xg - anchors(i,1)).^2 + (Yg - anchors(i,2)).^2);
    Vg = Vg + 0.5 * k * (dg - L0).^2;
end
% Mask outside ring
outside  = (Xg.^2 + Yg.^2) > R^2;
Vg(outside) = NaN;

contourf(Xg, Yg, Vg, 40, 'LineStyle','none');
colormap(hot); colorbar;
hold on;
plot(anchors(:,1), anchors(:,2), 'co', 'MarkerFaceColor','c', 'MarkerSize',8);
plot(x(1:seg:end), y(1:seg:end), 'b.', 'MarkerSize', 1);
title('Potential energy landscape + trajectory','FontSize',12);
axis equal; xlabel('x'); ylabel('y');

%  equations of motion
function ds = eom(~, s, m, k, L0, c, anchors)
    px = s(1); py = s(2);
    vx = s(3); vy = s(4);

    Fx = 0; Fy = 0;
    for i = 1:size(anchors,1)
        ax  = anchors(i,1);
        ay  = anchors(i,2);
        dx  = ax - px;
        dy  = ay - py;
        d   = sqrt(dx^2 + dy^2);
        if d > 1e-10                     
            F   = k * (d - L0) / d;    
            Fx  = Fx + F * dx;
            Fy  = Fy + F * dy;
        end
    end

    % Damping
    Fx = Fx - c * vx;
    Fy = Fy - c * vy;

    ds = [vx; vy; Fx/m; Fy/m];
end
