function MESH_ep = generate_mesh_SillTip_11layers(mesharg)

% Part of Open-GeoNabla, copyright GPLv3, 2018
% https://github.com/albansouche/Open-GeoNabla/
% Physics of Geological Processes (PGP) , The NJORD Centre, Dept of Geosciences, University of Oslo
% Author: Alban Souche

% 2D Mesh generation build upon Mutils-0.4-2 (Triangle mesh generator)


%%%%% FULL DOMAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_max		=  mesharg.width_mod/2;
x_min		= -x_max;
y_max		= mesharg.depth_top_layer;
y_min		= -(mesharg.tck_l1+mesharg.tck_l2+mesharg.tck_l3+mesharg.tck_l4+mesharg.tck_l5+mesharg.tck_sill+mesharg.tck_l7+mesharg.tck_l8+mesharg.tck_l9+mesharg.tck_l10+mesharg.tck_l11) + y_max;

%%% mesharg.seg_dx_max
mesharg.seg_dx_max = sqrt(2*mesharg.max_tri_area);


%%%% ELASTO-PLASTIC DOMAIN  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pointlist_ep   = [];
segmentlist_ep = [];
regionlist_ep  = [];

%BOX
x1 = x_min                   ; y1 = y_min;
x2 = x_max                   ; y2 = y_min;
x3 = x_max                   ; y3 = y_min+mesharg.tck_l11;
x4 = x_max                   ; y4 = y3+mesharg.tck_l10;
x5 = x_max                   ; y5 = y4+mesharg.tck_l9;
x6 = x_max                   ; y6 = y5+mesharg.tck_l8;
x7 = x_max                   ; y7 = y6+mesharg.tck_l7;
x8 = x_max                   ; y8 = y7+mesharg.tck_sill;
x9 = x_max                   ; y9 = y8+mesharg.tck_l5;
x10 = x_max                  ; y10 = y9+mesharg.tck_l4;
x11 = x_max                  ; y11 = y10+mesharg.tck_l3;
x12 = x_max                  ; y12 = y11+mesharg.tck_l2;
x13 = x_max                  ; y13 = y_max;
x14 = x_min                  ; y14 = y13;
x15 = x_min                  ; y15 = y12;
x16 = x_min                  ; y16 = y11;
x17 = x_min                  ; y17 = y10;
x18 = x_min                  ; y18 = y9;
x19 = x_min                  ; y19 = y8;
x20 = x_min                  ; y20 = y7;
x21 = x_min                  ; y21 = y6;
x22 = x_min                  ; y22 = y5;
x23 = x_min                  ; y23 = y4;
x24 = x_min                  ; y24 = y3;

x25 = x_min+mesharg.lth_sill ; y25 = y7;   % along the sill
x26 = x25                    ; y26 = y8;   % along the sill 

tmp = 3*mesharg.tck_sill + 1; %min(mesharg.lth_sill + (mesharg.width_mod-mesharg.lth_sill)/3, mesharg.lth_sill + 5*mesharg.tck_sill ); min(mesharg.lth_sill + (mesharg.width_mod-mesharg.lth_sill)/3, mesharg.lth_sill + 5*mesharg.tck_sill );
x27 = x_min +tmp              ; y27 = y11; % x limit of refinement
x28 = x_min +tmp              ; y28 = y10; % x limit of refinement
x29 = x_min +tmp              ; y29 = y9; % x limit of refinement
x30 = x_min +tmp              ; y30 = y8;  % x limit of refinement
x31 = x_min +tmp              ; y31 = y7;  % x limit of refinement
x32 = x_min +tmp              ; y32 = y6;  % x limit of refinement
x33 = x_min +tmp              ; y33 = y5;  % x limit of refinement 
x34 = x_min +tmp              ; y34 = y4;  % x limit of refinement 


BOX =     [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18 x19 x20 x21 x22 x23 x24 x25 x26 x27 x28 x29 x30 x31 x32 x33 x34 ;...
          y1 y2 y3 y4 y5 y6 y7 y8 y9 y10 y11 y12 y13 y14 y15 y16 y17 y18 y19 y20 y21 y22 y23 y24 y25 y26 y27 y28 y29 y30 y31 y32 y33  y34 ];
BOX_p  = [x13-1e-6  x13-1e-6  x13-1e-6  x13-1e-6  x13-1e-6  x13-1e-6 x13-1e-6  x13-1e-6  x13-1e-6  x13-1e-6 x13-1e-6 x19+1e-6  x27-1e-6 x28-1e-6 x29-1e-6  x32-1e-6 x33-1e-6;...
          y13-1e-6  y12-1e-6  y11-1e-6  y10-1e-6  y9-1e-6   y8-1e-6  y7-1e-6   y6-1e-6   y5-1e-6   y4-1e-6  y3-1e-6  y19-1e-6  y27-1e-6 y28-1e-6 y29-1e-6  y32-1e-6 y33-1e-6;...
                 1         2         3         4        5         6        7         8         9        10       11        12        13       14       15        16       17];

pointlist_ep  = [pointlist_ep   BOX];
regionlist_ep = [regionlist_ep  BOX_p];

% CREATE CRACK TIP AND INTRUSION
bar_arcang  = 180; % arc angle of the bar section on the disk
bar_arcang  = bar_arcang*pi/180;
theta       = linspace(pi-(bar_arcang/2) ,-pi+(bar_arcang/2) ,mesharg.no_pts_tip);
xx_incl     = (mesharg.tck_sill/2) * cos(theta);
yy_incl     = (mesharg.tck_sill/2) * sin(theta);
dx_incl     = sqrt( (xx_incl(2:end)-xx_incl(1:end-1)).^2 + (yy_incl(2:end)-yy_incl(1:end-1)).^2 );
indx_split = find(dx_incl>mesharg.seg_dx_max(6));
while ~isempty(indx_split)
    %update nodes
    new_x_nodes  = (xx_incl(2:end)+xx_incl(1:end-1)) /2;
    new_y_nodes  = (yy_incl(2:end)+yy_incl(1:end-1)) /2;
    xx_incl      = [ xx_incl; [new_x_nodes 0]];
    xx_incl      = xx_incl(:)';
    xx_incl(end) = [];
    yy_incl      = [ yy_incl; [new_y_nodes 0]];
    yy_incl      = yy_incl(:)';
    yy_incl(end) = [];
    %reevaluate the distance between segments
    dx_incl = sqrt( (xx_incl(2:end)-xx_incl(1:end-1)).^2 + (yy_incl(2:end)-yy_incl(1:end-1)).^2 );
    indx_split = find(dx_incl>mesharg.seg_dx_max(6));  
end
xx_incl([1 end]) = [];
yy_incl([1 end]) = [];
xx_incl        = xx_incl + x25;
yy_incl        = yy_incl + (y25+y26)/2;


% seg6 (13-14) (tip right)
TOADD    = [xx_incl; yy_incl];
no_pts   = size(TOADD,2);
no_pts_o   = size(pointlist_ep,2);
TOADD_s  = [26;no_pts_o+1;mesharg.seg_id(6)];
TOADD_s  = [TOADD_s [no_pts_o+(1:no_pts);(1:no_pts)+no_pts_o+1; mesharg.seg_id(6)*ones(1,no_pts)] ];
TOADD_s(2,end)  = 25;
pointlist_ep    = [pointlist_ep   TOADD];
segmentlist_ep  = [segmentlist_ep TOADD_s];

% seg1 bottom
TOADD_s  = [1   2 mesharg.seg_id(1) ; ...
            2   3 mesharg.seg_id(2) ; ...
            3   4 mesharg.seg_id(2) ; ...
            4   5 mesharg.seg_id(2) ; ...
            5   6 mesharg.seg_id(2) ; ...
            6   7 mesharg.seg_id(2) ; ...
            7   8 mesharg.seg_id(2) ; ...
            8   9 mesharg.seg_id(2) ; ...    
            9  10 mesharg.seg_id(2) ; ...
            10 11 mesharg.seg_id(2) ; ...
            11 12 mesharg.seg_id(2) ; ...
            12 13 mesharg.seg_id(2) ; ...
            13 14 mesharg.seg_id(3) ; ...
            14 15 mesharg.seg_id(4) ; ...
            15 16 mesharg.seg_id(4) ; ...
            16 17 mesharg.seg_id(4) ; ...
            17 18 mesharg.seg_id(4) ; ...
            18 19 mesharg.seg_id(4) ; ...
            19 20 mesharg.seg_id(5) ; ...
            20 21 mesharg.seg_id(4) ; ...
            21 22 mesharg.seg_id(4) ; ...
            22 23 mesharg.seg_id(4) ; ...
            23 24 mesharg.seg_id(4) ; ...
            24  1 mesharg.seg_id(4) ; ...
            26 19 mesharg.seg_id(6) ; ...
            25 20 mesharg.seg_id(6) ; ...
            24 3   11 ; ...
            23 34  11 ; ...
            34 4   11 ; ...
            22 33  11 ; ...
            33 5   11 ; ...
            21 32  11 ; ...
            32 6   11 ; ...
            31 7   11 ; ...
            30 8   11 ; ...
            18 29  11 ; ...
            29 9   11 ; ...
            17 28  11 ; ...
            28 10  11 ; ...
            16 27  11 ; ...
            27 11  11 ; ...
            15 12  11 ; ...
            34 33  11 ; ...
            33 32  11 ; ...
            32 29  11 ; ...
            29 28  11 ; ...
            28 27  11];    

segmentlist_ep  = [segmentlist_ep TOADD_s'];


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
opts.max_tri_area     = [mesharg.max_tri_area mesharg.max_tri_area(end)*ones(1,5)] ;     % maximum triangle area

opts.ignore_holes     = 0;     %
opts.exact_arithmetic = 1;     % 0: do not use Triangle's exact arithmetic
opts.zero_numbering   = 0;     % 1: use zero-based index numbering

%%% OPTIONS TETGEN
switch mesharg.type_el_ep 
    case 'tri3'
        opts.other_options   = 'pqAa';
    case 'tri7'
        opts.other_options   = 'po2qAa';
end

% Create triangle input structure
tristr.points         = pointlist_ep;
tristr.segments       = uint32(segmentlist_ep(1:2,:));  % note segments have to be uint32
tristr.segmentmarkers = uint32(segmentlist_ep(3,:));
tristr.regions        = [regionlist_ep; opts.max_tri_area ];

% figure(1), clf
% plot(tristr.points(1,:),tristr.points(2,:),'rd')
% hold on
% for i=1:length(tristr.segments)
% plot(tristr.points(1,tristr.segments(:,i)),tristr.points(2,tristr.segments(:,i)),'-k')
% drawnow
% end

% Generate the mesh using triangle
MESH      = mtriangle(opts, tristr);
MESH.opts = opts;

% Rename element number 13 14 15 16 17 18 19 to 3 4 5 6 7 8 9
MESH.elem_markers(MESH.elem_markers==13) = 3;
MESH.elem_markers(MESH.elem_markers==14) = 4;
MESH.elem_markers(MESH.elem_markers==15) = 6;
MESH.elem_markers(MESH.elem_markers==16) = 8;
MESH.elem_markers(MESH.elem_markers==17) = 9;

% Special treatment for region 5 and 7 
indx_ele = find(MESH.elem_markers==6);
indx = MESH.NODES(2,MESH.ELEMS(7,indx_ele)) > (mesharg.tck_sill/2);
MESH.elem_markers(indx_ele(indx)) = 5;
indx = MESH.NODES(2,MESH.ELEMS(7,indx_ele)) < (-mesharg.tck_sill/2);
MESH.elem_markers(indx_ele(indx)) = 7;


% Delete inner segmentation (segment nb 11)
indx_seg11 = MESH.segment_markers==11;
MESH.segment_markers(indx_seg11) = [];
MESH.SEGMENTS(:,indx_seg11)      = [];
% Additional field
% MESH.sill_interface_nodes = [];

toplot = repmat(MESH.elem_markers,3,1);
trisurf(reshape(1:3*size(MESH.ELEMS,2),3, size(MESH.ELEMS,2))', MESH.NODES(1,MESH.ELEMS(1:3,:)), MESH.NODES(2,MESH.ELEMS(1:3,:)), zeros(size(MESH.NODES(1,MESH.ELEMS(1:3,:)))), toplot);
title('mechanical mesh'); view(2); axis image; axis off
drawnow  
   
% SPLITTING MESH INTO ElatoPlastic AND Viscous DOMAINS

%%% FIND NEIGHBORGS IN MESH
if exist('MESH.NEIGHBORS','var')==0 
    tmp_elems = double(MESH.ELEMS(1:3,:));
    tmp_nodes = MESH.NODES(:,1:max(tmp_elems(:)));
    trep = triangulation(tmp_elems', tmp_nodes');
    MESH.NEIGHBORS = trep.neighbors()';
    MESH.NEIGHBORS(isnan(MESH.NEIGHBORS)) = 0;
    MESH.NEIGHBORS = uint32(MESH.NEIGHBORS);
end

%% MESHv for visous domain %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nb_fluid_dom = 12;
%
%

%% MESH_ep for elastoplastic domain %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MESH_ep.ELEMS = MESH.ELEMS(:,MESH.elem_markers<nb_fluid_dom);
[uniq_nodes, ~, perm_nodes] = unique(MESH_ep.ELEMS);
iperm_nodes                 = zeros(1,max(MESH_ep.ELEMS(:)));
iperm_nodes(MESH_ep.ELEMS)  = perm_nodes;
MESH_ep.NODES               = MESH.NODES(:,uniq_nodes);
MESH_ep.ELEMS               = uint32( reshape(perm_nodes, size(MESH_ep.ELEMS)) );
tmp_seg                     = ~(MESH.segment_markers==mesharg.seg_id(5));
MESH_ep.SEGMENTS            = MESH.SEGMENTS(:,tmp_seg);
MESH_ep.SEGMENTS            = iperm_nodes(MESH_ep.SEGMENTS);
MESH_ep.segment_markers     = MESH.segment_markers(tmp_seg);
MESH_ep.elem_markers        = MESH.elem_markers(MESH.elem_markers<nb_fluid_dom);

% update node markers based on facet numbering
indx_bt = find( MESH_ep.segment_markers==mesharg.seg_id(1) | MESH_ep.segment_markers==mesharg.seg_id(3) ); % index boundary facet (bottom/top)
indx_2  = find( MESH_ep.segment_markers==mesharg.seg_id(2) ); % index boundary right
indx_4  = find( MESH_ep.segment_markers==mesharg.seg_id(4) ); % index boundary left
indx_6  = find( MESH_ep.segment_markers==mesharg.seg_id(6) ); % intrusion body

MESH_ep.node_markers = zeros(1,size(MESH_ep.NODES,2));
nnps                 = 2; % number of segments per nodes
MESH_ep.node_markers(MESH_ep.SEGMENTS(:,indx_bt)) = repmat(MESH_ep.segment_markers(indx_bt),nnps,1);
MESH_ep.node_markers(MESH_ep.SEGMENTS(:,indx_2))  = repmat(MESH_ep.segment_markers(indx_2),nnps,1);
MESH_ep.node_markers(MESH_ep.SEGMENTS(:,indx_6))  = repmat(MESH_ep.segment_markers(indx_6),nnps,1);  % index boundary left
MESH_ep.node_markers(MESH_ep.SEGMENTS(:,indx_4))  = repmat(MESH_ep.segment_markers(indx_4),nnps,1);  % index boundary left

%  paraview_write('tet4',MESH_ep.NODES, MESH_ep.ELEMS, MESH_ep.elem_markers, MESH_ep.node_markers,'mesh_ep_ini', 0)

%%% FIND NEIGHBORGS IN MESH_ep
if exist('MESH_ep.NEIGHBORS','var')==0 
    tmp_elems = double(MESH_ep.ELEMS(1:3,:));
    tmp_nodes = MESH_ep.NODES(:,1:max(tmp_elems(:)));
    trep = triangulation(tmp_elems', tmp_nodes');
    MESH_ep.NEIGHBORS = trep.neighbors()';
    MESH_ep.NEIGHBORS(isnan(MESH_ep.NEIGHBORS)) = 0;
    MESH_ep.NEIGHBORS = uint32(MESH_ep.NEIGHBORS);
end

%% REODERING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spopts.symmetric = 0;
spopts.ndof      = 1;
K = sparse_create(MESH_ep.ELEMS,1, spopts);
% [perm,iperm] = mrcm(K); % MRCM reordering of NODES
perm = uint32(amd(K));               % AMD reordering of NODES
iperm = zeros(size(perm), 'uint32'); % AMD reordering of NODES
iperm(perm)  = 1:length(perm);       % AMD reordering of NODES

MESH_ep.NODES = MESH_ep.NODES(:,perm);
MESH_ep.node_markers = MESH_ep.node_markers(:,perm);
MESH_ep.SEGMENTS = reshape(iperm(MESH_ep.SEGMENTS), size(MESH_ep.SEGMENTS));
MESH_ep.ELEMS = reshape(iperm(MESH_ep.ELEMS), size(MESH_ep.ELEMS));

% Geometric Reordering of ELEMS
[~, permel]  = sort(min(MESH_ep.ELEMS)); 
% [~, permel]  = sort(mean(MESH_ep.ELEMS)); 
permel       = uint32(permel);
MESH_ep.ELEMS  = MESH_ep.ELEMS(:,permel);
MESH_ep.elem_markers = MESH_ep.elem_markers(permel);
if isfield(MESH_ep, 'NEIGHBORS')
    % first permute the elements
    MESH_ep.NEIGHBORS = MESH_ep.NEIGHBORS(:, permel);
    % now the neighbor information for every element
    noneighbor = (MESH_ep.NEIGHBORS==0);
    MESH_ep.NEIGHBORS(noneighbor) = 1;
    ipermel(permel)= uint32(1:length(permel));
    MESH_ep.NEIGHBORS = ipermel(MESH_ep.NEIGHBORS);
    MESH_ep.NEIGHBORS(noneighbor) = 0;
end

% Calculate the normals to the external segments
MESH_ep = Normal_calc_2D(MESH_ep);



