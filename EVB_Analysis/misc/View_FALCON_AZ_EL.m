close all; clear all; clc;
% --- User Parameters ---
az_deg = -12;      % Azimuth angle in degrees (rotation about Z axis)
el_deg = 24;      % Elevation angle in degrees (rotation about Y axis)
stl_file = 'falcon_model.stl';  % STL file name

% --- Load STL ---
model = stlread(stl_file);  % returns a triangulation object

% Get faces and vertices
F = model.ConnectivityList;
V = model.Points;

% --- Convert angles to radians ---
az = deg2rad(az_deg);
el = deg2rad(el_deg);

% --- Define rotation matrices ---
Rz = [cos(az) -sin(az) 0;
      sin(az)  cos(az) 0;
         0       0     1];  % Azimuth: about Z

Ry = [cos(el)  0 sin(el);
         0     1    0;
     -sin(el) 0 cos(el)];  % Elevation: about Y

% --- Apply combined rotation ---
R = Rz * Ry;         % First rotate about Y (elevation), then Z (azimuth)
V_rot = (R * V')';   % Apply rotation to all vertices

% --- Plotting ---
figure('Name', 'FALCON Rotated Views', 'Position', [100 100 1000 800]);
sgtitle(sprintf('Rotated Model: Az = %d°, El = %d°', az_deg, el_deg));

% View 1: XY (Top-Down)
subplot(2,2,1); hold on;
patch('Faces', F, 'Vertices', V_rot, ...
      'FaceColor', [0.6 0.8 1], 'EdgeColor', 'k');
view(0, 90);  % XY plane
axis equal; title('Top View');
xlabel('X'); ylabel('Y'); zlabel('Z');
grid on;

% View 2: XZ (Side View 1)
subplot(2,2,2); hold on;
patch('Faces', F, 'Vertices', V_rot, ...
      'FaceColor', [0.6 0.8 1], 'EdgeColor', 'k');
view(0, 0);   % Looking along +Y
axis equal; title('Side View');
xlabel('X'); ylabel('Y'); zlabel('Z');
grid on;

% View 3: XZ (Side View 2, opposite side)
subplot(2,2,3); hold on;
patch('Faces', F, 'Vertices', V_rot, ...
      'FaceColor', [0.6 0.8 1], 'EdgeColor', 'k');
view(180, 0); % Looking along -Y
axis equal; title('Side View');
xlabel('X'); ylabel('Y'); zlabel('Z');
grid on;

% View 4: YZ (Front View)
subplot(2,2,4); hold on;
patch('Faces', F, 'Vertices', V_rot, ...
      'FaceColor', [0.6 0.8 1], 'EdgeColor', 'k');
view(90, 0);  % Looking along +X
axis equal; title('Front View');
xlabel('X'); ylabel('Y'); zlabel('Z');
grid on;

% --- Lighting for all ---
for i = 1:4
    subplot(2,2,i);
    camlight headlight; lighting gouraud;
end
