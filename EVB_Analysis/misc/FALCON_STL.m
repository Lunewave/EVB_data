% Dimensions
plate_size = 508;        % mm
plate_thickness = 5;     % mm
prong_height = 30;       % mm
special_prong_height = 60; % Taller special prong
standard_radius = 52;    % mm
special_radius  = 72;    % mm
num_prongs = 5;
cyl_radius = 1;          % mm

% Centering offset
offset_xy = -plate_size / 2;
offset_z  = -plate_thickness / 2;

figure; hold on;

% Define and shift plate vertices
plate_vertices = [
    0 0 0;
    plate_size 0 0;
    plate_size plate_size 0;
    0 plate_size 0;
    0 0 plate_thickness;
    plate_size 0 plate_thickness;
    plate_size plate_size plate_thickness;
    0 plate_size plate_thickness;
];
plate_vertices(:,1:2) = plate_vertices(:,1:2) + offset_xy;
plate_vertices(:,3)   = plate_vertices(:,3)   + offset_z;

% Flip X coordinates to mirror figure about YZ plane
plate_vertices(:,1) = -plate_vertices(:,1);

plate_faces = [
    1 2 3 4;
    5 6 7 8;
    1 2 6 5;
    2 3 7 6;
    3 4 8 7;
    4 1 5 8;
];

patch('Vertices', plate_vertices, 'Faces', plate_faces, ...
      'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'k');

% Plate center (now origin)
cx = 0;
cy = 0;
cz = 0 + plate_thickness / 2;

% Prong resolution
n_circ = 30;
[XC, YC, ZC] = cylinder(cyl_radius, n_circ);
ZC_std = ZC * prong_height;
ZC_special = ZC * special_prong_height;

% --- Standard Prongs 1 to 5 ---
for k = 0:num_prongs-1
    theta = deg2rad(k * 72);
    % Flip X direction by subtracting cosine term
    px = cx - standard_radius * cos(theta);
    py = cy + standard_radius * sin(theta);

    X = XC + px;
    Y = YC + py;
    Z = ZC_std + cz;

    surf(X, Y, Z, 'FaceColor', [0.2 0.5 1], 'EdgeColor', 'none');

    % Number label
    text(px, py, cz + prong_height + 5, ...
         sprintf('%d', k+1), ...
         'FontSize', 14, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'center', 'Color', 'k');
end

% --- Special 6th Prong ---
% Parameters for L-prong
L_base_radius = 72;   % mm, distance from center
L_leg_height = 90;    % vertical leg height
L_leg_thickness = 2;  % cross-section of beam (square)
L_leg_length = 72;    % length of horizontal leg pointing inward

% Flip X direction here too
base_x = cx - L_base_radius;  % flipped
base_y = cy;
base_z = cz;

% --- Vertical leg of L --- (rectangular prism)
vert_verts = [
    0  -L_leg_thickness/2 0;
    L_leg_thickness  -L_leg_thickness/2 0;
    L_leg_thickness  L_leg_thickness/2 0;
    0  L_leg_thickness/2 0;
    0  -L_leg_thickness/2 L_leg_height;
    L_leg_thickness  -L_leg_thickness/2 L_leg_height;
    L_leg_thickness  L_leg_thickness/2 L_leg_height;
    0  L_leg_thickness/2 L_leg_height;
];

% Shift vertical leg to base position
vert_verts(:,1) = vert_verts(:,1) + base_x;
vert_verts(:,2) = vert_verts(:,2) + base_y;
vert_verts(:,3) = vert_verts(:,3) + base_z;

% Faces of a box (quads)
box_faces = [
    1 2 3 4;  % bottom
    5 6 7 8;  % top
    1 2 6 5;
    2 3 7 6;
    3 4 8 7;
    4 1 5 8;
];

patch('Vertices', vert_verts, 'Faces', box_faces, ...
      'FaceColor', [1 0.4 0.2], 'EdgeColor', 'k');

% --- Horizontal leg of L --- (rectangular prism pointing inward, along -X)
horz_verts = [
    L_leg_length  -L_leg_thickness/2 0;
    0             -L_leg_thickness/2 0;
    0             L_leg_thickness/2 0;
    L_leg_length  L_leg_thickness/2 0;
    L_leg_length  -L_leg_thickness/2 L_leg_thickness;
    0             -L_leg_thickness/2 L_leg_thickness;
    0             L_leg_thickness/2 L_leg_thickness;
    L_leg_length  L_leg_thickness/2 L_leg_thickness;
];

% Shift horizontal leg to top of vertical leg base (same base_x,y,z)
horz_verts(:,1) = horz_verts(:,1) + base_x;
horz_verts(:,2) = horz_verts(:,2) + base_y;
horz_verts(:,3) = horz_verts(:,3) + base_z + L_leg_height - L_leg_thickness;

patch('Vertices', horz_verts, 'Faces', box_faces, ...
      'FaceColor', [1 0.4 0.2], 'EdgeColor', 'k');

% Add label for prong 6
text(base_x, base_y, base_z + L_leg_height + 5, ...
     '6', 'FontSize', 14, 'FontWeight', 'bold', ...
     'HorizontalAlignment', 'center', 'Color', 'k');

% Small rectangular prism dimensions
small_thickness = 40;  % Z direction height
small_height = 2;      % Y direction thickness
small_length = 15;     % X direction length

% Define vertices for the small prism
small_verts = [
    0                -small_height/2    0;
    small_length      -small_height/2    0;
    small_length       small_height/2    0;
    0                 small_height/2    0;
    0                -small_height/2    small_thickness;
    small_length      -small_height/2    small_thickness;
    small_length       small_height/2    small_thickness;
    0                 small_height/2    small_thickness;
];

% Position the small prism at the end of the horizontal leg
small_verts(:,1) = small_verts(:,1) + base_x + L_leg_length;
small_verts(:,2) = small_verts(:,2) + base_y;
small_verts(:,3) = small_verts(:,3) + base_z + L_leg_height - L_leg_thickness - small_thickness/2;

% Plot the small rectangular prism
patch('Vertices', small_verts, 'Faces', box_faces, ...
      'FaceColor', [0.8 0.2 0.2], 'EdgeColor', 'k');

% --- Helper function to convert quads to triangles ---
function triFaces = quads2triangles(quadFaces)
    triFaces = zeros(size(quadFaces,1)*2, 3);
    for i = 1:size(quadFaces,1)
        q = quadFaces(i,:);
        triFaces(2*i-1,:) = q([1 2 3]);
        triFaces(2*i,:)   = q([1 3 4]);
    end
end

% --- Prepare STL export ---

% Convert quad faces to triangles
plate_faces_tri = quads2triangles(plate_faces);
box_faces_tri = quads2triangles(box_faces);

all_V = [];
all_F = [];
f_offset = 0;

% Plate
all_V = [all_V; plate_vertices];
all_F = [all_F; plate_faces_tri + f_offset];
f_offset = size(all_V,1);

% Standard prongs
for k = 0:num_prongs-1
    theta = deg2rad(k * 72);
    px = cx - standard_radius * cos(theta);
    py = cy + standard_radius * sin(theta);

    X = XC + px;
    Y = YC + py;
    Z = ZC_std + cz;

    [F,V] = surf2patch(X, Y, Z, 'triangles');

    all_V = [all_V; V];
    all_F = [all_F; F + f_offset];
    f_offset = size(all_V,1);
end

% Vertical leg
all_V = [all_V; vert_verts];
all_F = [all_F; box_faces_tri + f_offset];
f_offset = size(all_V,1);

% Horizontal leg
all_V = [all_V; horz_verts];
all_F = [all_F; box_faces_tri + f_offset];
f_offset = size(all_V,1);

% Small prism
all_V = [all_V; small_verts];
all_F = [all_F; box_faces_tri + f_offset];
f_offset = size(all_V,1);

% Write STL file
% Create triangulation object and export to STL
TR = triangulation(all_F, all_V);
stlwrite(TR, 'falcon_model.stl');


% Final view settings
axis equal;
xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');
view(30, 30);
camlight; lighting gouraud;
title('FALCON');
