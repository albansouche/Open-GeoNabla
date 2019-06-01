function MESH_ep = generate_mesh_Gerbault2012(mesharg)

% Part of Open-GeoNabla, copyright GPLv3, 2018
% https://github.com/albansouche/Open-GeoNabla/
% Physics of Geological Processes (PGP) , The NJORD Centre, Dept of Geosciences, University of Oslo
% Author: Alban Souche

% 2D Mesh generation build upon Mutils-0.4-2 (Triangle mesh generator)


%% FULL DOMAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_min = 0;
x_max = mesharg.Lx;
y_min = - mesharg.Ly;
y_max = 0;
y_mag_mid = - mesharg.Mag_depth;
y_mag_max = - mesharg.Mag_depth + mesharg.Mag_rad;
y_mag_min = - mesharg.Mag_depth - mesharg.Mag_rad;
x_ref = mesharg.x_ref;
y_ref = - mesharg.y_ref;

%%%% ELASTO-PLASTIC DOMAIN  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pointlist_ep   = [];
segmentlist_ep = [];
regionlist_ep  = [];

%% BOX
BOX   = [x_min x_max x_max x_ref x_min x_min     x_min     x_min x_ref  ;...
         y_min y_min y_max y_max y_max y_mag_max y_mag_min y_ref y_ref ];
BOX_s = [1 2 3 4 5 7 8 8 4;...
         2 3 4 5 6 8 1 9 9;...
         mesharg.seg_id(1) mesharg.seg_id(2) mesharg.seg_id(3) mesharg.seg_id(3) mesharg.seg_id(4) mesharg.seg_id(4) mesharg.seg_id(4) 11 11];

BOX_p = [x_min+1e-3 x_min+1e-3;y_min+1e-3 y_ref+1e-3; 1 2];

pointlist_ep    = [pointlist_ep   BOX];
segmentlist_ep  = [segmentlist_ep BOX_s];
regionlist_ep   = [regionlist_ep  BOX_p];
no_pts_o_ep     = size(pointlist_ep,2);


%% CREATE CRACK TIP AND INTRUSION
theta    = linspace(pi/2,-pi/2,mesharg.Mag_nb_pts);
xx_incl = mesharg.Mag_rad * cos(theta);
yy_incl = mesharg.Mag_rad * sin(theta) + y_mag_mid;
xx_incl([1 end]) = [];
yy_incl([1 end]) = [];

% hold on
% plot(xx_incl, yy_incl,'rd')
TOADD    = [xx_incl; yy_incl];
no_pts   = size(TOADD,2);
TOADD_s  = [6;no_pts_o_ep+1;mesharg.seg_id(6)];  
TOADD_s  = [TOADD_s [no_pts_o_ep+1:no_pts+no_pts_o_ep;(no_pts_o_ep+1:no_pts+no_pts_o_ep)+1; mesharg.seg_id(6)*ones(1,no_pts)] ];
TOADD_s(2,end)  = 7;

pointlist_ep    = [pointlist_ep   TOADD];
segmentlist_ep  = [segmentlist_ep TOADD_s];


%%%%%%%%%%% mtriangle from mutils
% Set triangle options
opts = [];
opts.element_type     = mesharg.type_el_ep;% use 'tri3', 'tri6', or 'tri7'
opts.triangulate_poly = 1;     % Triangulates a Planar Straight Line Graph
opts.gen_edges        = 0;     % 1: return edge information
opts.gen_neighbors    = 1;     % 1: return element neighbor information
opts.gen_elmarkers    = 1;     % 1: return element markers
opts.gen_boundaries   = 1;     % 1: return node markers
opts.min_angle        = 30;    % minimum triangle angle
opts.max_tri_area     = mesharg.max_tri_area;     % maximum triangle area
opts.ignore_holes     = 1;     %
opts.exact_arithmetic = 1;     % 0: do not use Triangle's exact arithmetic
opts.zero_numbering   = 0;     % 1: use zero-based index numbering
opts.other_options    = 'pa';    % other triangle options

% Create triangle input structure
tristr.points         = pointlist_ep;
tristr.segments       = uint32(segmentlist_ep(1:2,:));  % note segments have to be uint32
tristr.segmentmarkers = uint32(segmentlist_ep(3,:));
tristr.regions        = [regionlist_ep; mesharg.max_tri_area ];

% for i=1:length(tristr.segments )
% plot(tristr.points(1,tristr.segments(:,i)), tristr.points(2,tristr.segments(:,i)),'-gd')
% drawnow
% pause
% end

% Generate the mesh using triangle
MESH_ep      = mtriangle(opts, tristr);
MESH_ep.opts = opts;

% Additional manupulation
MESH_ep.sill_interface_nodes = [];

MESH_ep.elem_markers = 1 + 0*MESH_ep.elem_markers; 

MESH_ep.SEGMENTS(:,MESH_ep.segment_markers==11) = [];
MESH_ep.segment_markers(:,MESH_ep.segment_markers==11) = [];


