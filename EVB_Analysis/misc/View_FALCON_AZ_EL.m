close all; clear all; clc;

% --- User Parameters ---
az_deg = 36;      % Azimuth angle in degrees (rotation about Z axis)
el_deg = 0;       % Elevation angle in degrees (rotation about Y axis)
stl_file = 'falcon_model.stl';  % STL file name

% --- Load STL ---
model = stlread(stl_file);  % triangulation object
F = model.ConnectivityList;
V = model.Points;

% --- Convert angles to radians ---
az = deg2rad(az_deg);
el = deg2rad(el_deg);

% --- Define rotation matrices ---
Rz = [cos(az), -sin(az), 0;
      sin(az),  cos(az), 0;
         0,        0,    1];   % Azimuth (Z)
     
Ry = [cos(el),  0, sin(el);
         0,     1,   0;
     -sin(el), 0, cos(el)];    % Elevation (Y)

% --- Apply rotation: Rz * Ry
R = Rz * Ry;
V_rot = (R * V')';

% --- Create mirrored version (flip across Y and Z only, NOT X) ---
V_mirror = V;
V_mirror(:,2) = -V_mirror(:,2);  % mirror Y
V_mirror(:,3) = -V_mirror(:,3);  % mirror Z


% --- Define rotation matrices ---
Rz = [cos(az), -sin(az), 0;
      sin(az),  cos(az), 0;
         0,        0,    1];   % Azimuth (Z)
     
Ry = [cos(el),  0, sin(el);
         0,     1,   0;
     -sin(el), 0, cos(el)];    % Elevation (Y)

% --- Apply rotation: Rz * Ry
R = Rz * Ry;
V_mirror = (R * V_mirror')';




% --- Plotting ---
figure('Name', 'FALCON Rotated & Mirrored Views', 'Position', [100 100 1000 800]);
sgtitle(sprintf('Rotated & Mirrored Model: Az = %d°, El = %d°', az_deg, el_deg));

% View 1: XY (Top-Down)
subplot(2,2,1); hold on;
patch('Faces', F, 'Vertices', V_rot, ...
      'FaceColor', [0.6 0.8 1], 'EdgeColor', 'k');
patch('Faces', F, 'Vertices', V_mirror, ...
      'FaceColor', [1 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
% Fixed grid for vector field (adjust size to suit your STL)
[xg, yg, zg] = meshgrid(500:-100:-500, -100:100:100, -100:25:100);  % long X path from +500 to -500

% All arrows point in -X direction
u = -ones(size(xg));
v = zeros(size(yg));
w = zeros(size(zg));

% Plot arrows
quiver3(xg, yg, zg, u, v, w, 0.2, 'Color', 'r', 'LineWidth', 1.2, 'MaxHeadSize', 2);
view(az_deg, 0);  % Top-down view
axis equal; title('Mirror Image');
xlabel('X'); ylabel('Y'); zlabel('Z'); grid on;

% View 2: XZ (Side View 1)
subplot(2,2,2); hold on;
patch('Faces', F, 'Vertices', V_rot, ...
      'FaceColor', [0.6 0.8 1], 'EdgeColor', 'k');
view(0, 0);  % Looking along +Y
axis equal; title('Side View (XZ, +Y)');
xlabel('X'); ylabel('Y'); zlabel('Z'); grid on;

% View 3: XZ (Side View 2)
subplot(2,2,3); hold on;
patch('Faces', F, 'Vertices', V_rot, ...
      'FaceColor', [0.6 0.8 1], 'EdgeColor', 'k');
view(180, 0);  % Looking along -Y
axis equal; title('Side View (XZ, -Y)');
xlabel('X'); ylabel('Y'); zlabel('Z'); grid on;

% View 4: YZ (Front View)
subplot(2,2,4); hold on;
patch('Faces', F, 'Vertices', V_rot, ...
      'FaceColor', [0.6 0.8 1], 'EdgeColor', 'k');
view(90, 0);  % Looking along +X
axis equal; title('Front View (YZ)');
xlabel('X'); ylabel('Y'); zlabel('Z'); grid on;

% Add lighting to each subplot
for i = 1:4
    subplot(2,2,i);
    camlight headlight; lighting gouraud;
end

figure('Name', 'Selected Views Only (1 & 4)', 'Position', [1200 100 800 800]);
sgtitle(sprintf('Rotated & Mirrored Model: Az = %d°, El = %d°', az_deg, el_deg));


% --- Subplot 1 (Rotated Az/El View) ---
subplot(2,1,1); hold on;
patch('Faces', F, 'Vertices', V_rot, ...
      'FaceColor', [0.6 0.8 1], 'EdgeColor', 'k');
patch('Faces', F, 'Vertices', V_mirror, ...
      'FaceColor', [1 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
% Fixed grid for vector field (adjust size to suit your STL)
[xg, yg, zg] = meshgrid(500:-100:-500, -100:100:100, -100:25:100);  % long X path from +500 to -500

% All arrows point in -X direction
u = -ones(size(xg));
v = zeros(size(yg));
w = zeros(size(zg));

% Plot arrows
quiver3(xg, yg, zg, u, v, w, 0.2, 'Color', 'r', 'LineWidth', 1.2, 'MaxHeadSize', 2);
view(az_deg, 0);  % Top-down view
axis equal; title('Mirror Image');
xlabel('X'); ylabel('Y'); zlabel('Z'); grid on;

% --- Subplot 4 (YZ Front View) ---
subplot(2,1,2); hold on;
patch('Faces', F, 'Vertices', V_rot, ...
      'FaceColor', [0.6 0.8 1], 'EdgeColor', 'k');
view(90, 0);  % Looking along +X
axis equal; title('Front View (YZ)');
xlabel('X'); ylabel('Y'); zlabel('Z'); grid on;

% Add lighting
for i = 1:2
    subplot(2,1,i);
    camlight headlight; lighting gouraud;
end

