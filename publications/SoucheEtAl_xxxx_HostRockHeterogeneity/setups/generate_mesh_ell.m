function MESH = generate_mesh_ell(mesharg)

% Part of Open-GeoNabla, copyright GPLv3, 2018
% https://github.com/albansouche/Open-GeoNabla/
% Physics of Geological Processes (PGP) , The NJORD Centre, Dept of Geosciences, University of Oslo
% Author: Alban Souche

% 2D Mesh generation build upon Mutils-0.4-2 (Triangle mesh generator)


%%%%% FULL DOMAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_max		= mesharg.X_mod;
x_min		= 0;
y_max		= mesharg.Y_mod;
y_min		= 0;

%%%% ELASTO-PLASTIC DOMAIN  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pointlist_ep   = [];
segmentlist_ep = [];
regionlist_ep  = [];

%BOX
x1 = mesharg.ell_tck*mesharg.ell_ratio ; y1 = y_min;
x2 = mesharg.X_ref           ; y2 = y_min;
x3 = x_max                   ; y3 = y_min;
x4 = x_max                   ; y4 = y_max;
x5 = x_min                   ; y5 = y_max;
x6 = x_min                   ; y6 = mesharg.Y_ref;    %limit of refinement
x7 = x_min                   ; y7 = mesharg.ell_tck; %limit of refinement
x8 = mesharg.X_ref           ; y8 = mesharg.Y_ref;    %limit of refinement

BOX = [ x1 x2 x3 x4 x5 x6 x7 x8 ;...
        y1 y2 y3 y4 y5 y6 y7 y8 ];
BOX_p  = [ x1+1e-1  x4-1e-1 ;...
           y1+1e-1  y4-1e-1 ;...
                 2        1 ];

pointlist_ep  = [pointlist_ep   BOX];
regionlist_ep = [regionlist_ep  BOX_p];

% CREATE CRACK TIP AND INTRUSION
bar_arcang  = 90; % arc angle of the bar section on the disk
bar_arcang  = bar_arcang*pi/180;
theta       = linspace(bar_arcang ,0 ,mesharg.no_pts_tip);
xx_incl     = (mesharg.ell_tck*mesharg.ell_ratio) * cos(theta);
yy_incl     = mesharg.ell_tck * sin(theta);
xx_incl([1 end]) = [];
yy_incl([1 end]) = [];

% seg6 (intrusion segmentation)
TOADD    = [xx_incl; yy_incl];
no_pts   = size(TOADD,2);
no_pts_o   = size(pointlist_ep,2);
TOADD_s  = [7;no_pts_o+1;mesharg.seg_id(6)];
TOADD_s  = [TOADD_s [no_pts_o+(1:no_pts);(1:no_pts)+no_pts_o+1; mesharg.seg_id(6)*ones(1,no_pts)] ];
TOADD_s(2,end)  = 1;
pointlist_ep    = [pointlist_ep   TOADD];
segmentlist_ep  = [segmentlist_ep TOADD_s];

% seg1234 bottom
TOADD_s  = [1   2 mesharg.seg_id(1) ; ...
            2   3 mesharg.seg_id(1) ; ...
            3   4 mesharg.seg_id(2) ; ...
            4   5 mesharg.seg_id(3) ; ...
            5   6 mesharg.seg_id(4) ; ...
            6   7 mesharg.seg_id(4) ; ...
            6   8 11 ; ...
            8   2 11 ];

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
opts.max_tri_area     = mesharg.max_tri_area;     % maximum triangle area

opts.ignore_holes     = 1;     %
opts.exact_arithmetic = 1;     % 0: do not use Triangle's exact arithmetic
opts.zero_numbering   = 0;     % 1: use zero-based index numbering

%%% OPTIONS TETGEN
switch mesharg.type_el_ep 
    case 'tri3'
        opts.other_options   = 'pqAa';
    case 'tri7'
        opts.other_options   = 'po2qqAa';
end

% Create triangle input structure
tristr.points         = pointlist_ep;
tristr.segments       = uint32(segmentlist_ep(1:2,:));  % note segments have to be uint32
tristr.segmentmarkers = uint32(segmentlist_ep(3,:));
tristr.regions        = [regionlist_ep; [ opts.max_tri_area(2) opts.max_tri_area(1) ] ];

figure(1), clf
plot(tristr.points(1,:),tristr.points(2,:),'rd')
% hold on
% for i=1:length(tristr.segments)
% plot(tristr.points(1,tristr.segments(:,i)),tristr.points(2,tristr.segments(:,i)),'-k')
% drawnow
% end

% Generate the mesh using triangle
MESH      = mtriangle(opts, tristr);
MESH.opts = opts;

% Rename element number 9 10 11 to 3 4 5
MESH.elem_markers(MESH.elem_markers==2)  = 1;
% Delete inner segmentation (segment nb 11)
indx_seg11 = MESH.segment_markers==11;
MESH.segment_markers(indx_seg11) = [];
MESH.SEGMENTS(:,indx_seg11)      = [];
% Additional field
% MESH.sill_interface_nodes = [];

% figure(1),hold on
% toplot = repmat(MESH.elem_markers,3,1);
% trisurf(reshape(1:3*size(MESH.ELEMS,2),3, size(MESH.ELEMS,2))', MESH.NODES(1,MESH.ELEMS(1:3,:)), MESH.NODES(2,MESH.ELEMS(1:3,:)), zeros(size(MESH.NODES(1,MESH.ELEMS(1:3,:)))), toplot);
% title('mechanical mesh'); view(2); axis image; axis off
% drawnow  
   
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

% update node markers based on facet numbering
indx_bt = find( MESH.segment_markers==mesharg.seg_id(1) | MESH.segment_markers==mesharg.seg_id(3) ); % index boundary facet (bottom/top)
indx_2  = find( MESH.segment_markers==mesharg.seg_id(2) ); % index boundary right
indx_4  = find( MESH.segment_markers==mesharg.seg_id(4) ); % index boundary left
indx_6  = find( MESH.segment_markers==mesharg.seg_id(6) ); % intrusion body

nnps              = 2; % number of segments per nodes
MESH.node_markers = zeros(1,size(MESH.NODES,2));
MESH.node_markers(MESH.SEGMENTS(:,indx_bt)) = repmat(MESH.segment_markers(indx_bt),nnps,1);
MESH.node_markers(MESH.SEGMENTS(:,indx_2))  = repmat(MESH.segment_markers(indx_2),nnps,1);
MESH.node_markers(MESH.SEGMENTS(:,indx_4))  = repmat(MESH.segment_markers(indx_4),nnps,1);  % index boundary left
MESH.node_markers(MESH.SEGMENTS(:,indx_6))  = repmat(MESH.segment_markers(indx_6),nnps,1);  % index boundary left

%% REODERING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spopts.symmetric = 0;
spopts.ndof      = 1;
K = sparse_create(MESH.ELEMS,1, spopts);
% [perm,iperm] = mrcm(K); % MRCM reordering of NODES
perm = uint32(amd(K));               % AMD reordering of NODES
iperm = zeros(size(perm), 'uint32'); % AMD reordering of NODES
iperm(perm)  = 1:length(perm);       % AMD reordering of NODES

MESH.NODES = MESH.NODES(:,perm);
MESH.node_markers = MESH.node_markers(:,perm);
MESH.SEGMENTS = reshape(iperm(MESH.SEGMENTS), size(MESH.SEGMENTS));
MESH.ELEMS = reshape(iperm(MESH.ELEMS), size(MESH.ELEMS));

% Geometric Reordering of ELEMS
[~, permel]  = sort(min(MESH.ELEMS)); 
% [~, permel]  = sort(mean(MESH.ELEMS)); 
permel       = uint32(permel);
MESH.ELEMS  = MESH.ELEMS(:,permel);
MESH.elem_markers = MESH.elem_markers(permel);
if isfield(MESH, 'NEIGHBORS')
    % first permute the elements
    MESH.NEIGHBORS = MESH.NEIGHBORS(:, permel);
    % now the neighbor information for every element
    noneighbor = (MESH.NEIGHBORS==0);
    MESH.NEIGHBORS(noneighbor) = 1;
    ipermel(permel)= uint32(1:length(permel));
    MESH.NEIGHBORS = ipermel(MESH.NEIGHBORS);
    MESH.NEIGHBORS(noneighbor) = 0;
end

% Calculate the normals to the external segments
MESH = Normal_calc_2D(MESH);


