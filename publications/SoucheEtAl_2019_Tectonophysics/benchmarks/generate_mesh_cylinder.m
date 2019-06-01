function MESH = generate_mesh_cylinder(mesharg)

% Part of Open-GeoNabla, copyright GPLv3, 2018
% https://github.com/albansouche/Open-GeoNabla/
% Physics of Geological Processes (PGP) , The NJORD Centre, Dept of Geosciences, University of Oslo
% Author: Alban Souche

% 2D Mesh generation build upon Mutils-0.4-2 (Triangle mesh generator)


%%% Initiallization %%%%%%%%%%%%%
dx_av  = sqrt(2*mesharg.max_tri_area); % rought estimate 

%%% CREATE inner quarter disk %%%%%%%%%%%%
nb_pts_in = round( mesharg.angle_sym * mesharg.radius_in / dx_av );
theta_in = linspace(0 ,mesharg.angle_sym , nb_pts_in);
xx_in = mesharg.radius_in * cos(theta_in);
yy_in = mesharg.radius_in * sin(theta_in);

%%% CREATE external quarter disk %%%%%%%%%%%%
nb_pts_ext = round( mesharg.angle_sym * mesharg.radius_ext / dx_av );
theta_ext  = linspace(0 ,mesharg.angle_sym , nb_pts_ext);
xx_ext = mesharg.radius_ext * cos(theta_ext);
yy_ext = mesharg.radius_ext * sin(theta_ext);

%%% CREATE the upper segment %%%%%%%%%%
nb_pts_seg = round( (mesharg.radius_ext - mesharg.radius_in) / dx_av ); 
xx_seg_up = linspace(xx_in(end), xx_ext(end), nb_pts_seg);
yy_seg_up = linspace(yy_in(end), yy_ext(end), nb_pts_seg);

%%% CREATE the upper segment %%%%%%%%%%
xx_seg_low = linspace(xx_in(1), xx_ext(1), nb_pts_seg);
yy_seg_low = linspace(yy_in(1), yy_ext(1), nb_pts_seg);


%%% ASSEMBLING THE LISTS %%%%%%%%%%%
pointlist   = [xx_in xx_seg_up(2:end-1) fliplr(xx_ext) fliplr(xx_seg_low(2:end-1)) ;
               yy_in yy_seg_up(2:end-1) fliplr(yy_ext) fliplr(yy_seg_low(2:end-1)) ];
nb_nod      = length(pointlist);
segmentlist = [ 1:nb_nod ; [2:nb_nod 1]; [ ones(1,length(xx_in)-1) 2*ones(1,length(xx_seg_up)-1) 3*ones(1,length(xx_ext)-1)  4*ones(1,length(xx_seg_low)-1) ] ];
regionlist  = [mean(pointlist(1,:)); mean(pointlist(2,:)) ; 1];


%%%%%%%%%%% mtriangle from mutils
% Set triangle options
opts = [];
opts.element_type     = mesharg.type_el;% use 'tri3', 'tri6', or 'tri7'
opts.triangulate_poly = 1;     % Triangulates a Planar Straight Line Graph
opts.gen_edges        = 0;     % 1: return edge information
opts.gen_neighbors    = 1;     % 1: return element neighbor information
opts.gen_elmarkers    = 1;     % 1: return element markers
opts.gen_boundaries   = 1;     % 1: return node markers
opts.min_angle        = 30;    % minimum triangle angle
opts.max_tri_area     = mesharg.max_tri_area;     % maximum triangle area
opts.ignore_holes     = 0;     %
opts.exact_arithmetic = 1;     % 0: do not use Triangle's exact arithmetic
opts.zero_numbering   = 0;     % 1: use zero-based index numbering
opts.other_options    = '';    % other triangle options

 
% Create triangle input structure
tristr.points         = pointlist;
tristr.segments       = uint32(segmentlist(1:2,:));  % note segments have to be uint32
tristr.segmentmarkers = uint32(segmentlist(3,:));
tristr.regions        = [regionlist; opts.max_tri_area*ones(1,size(regionlist,2)) ];

% Generate the mesh using triangle
MESH      = mtriangle(opts, tristr);
MESH.opts = opts;


