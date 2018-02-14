function MESH = generate3Dmesh_ellipsoid(mesharg, SOLVER, modelname)

% Part of OpenGeoNabla, copyright GPLv3, 2018
% https://github.com/albansouche/OpenGeoNabla/
% Physics of Geological Processes (PGP) , The NJORD Centre, Dept of Geosciences, University of Oslo
% Author: Alban Souche


% 3D Mesh generation build upon Tetgen1.5.0 mesh generator


%%%%%%%%%%%%%% GENERATE PTS AND FACES OF DOMAIN BOX %%%%%%%%%%%%%%%%%%%%%%%

% Points
nb_point_box_edge = 16;
XYZ = zeros(nb_point_box_edge,3);

%main box
x_max = mesharg.width_mod/2;
x_min = -x_max;
y_max =  x_max;
y_min =  x_min;
z_max = 0;
z_min = -mesharg.depth_mod;

XYZ(1,1) = x_min;  XYZ(1,2) = y_max;  XYZ(1,3) = z_min;
XYZ(2,1) = x_max;  XYZ(2,2) = y_max;  XYZ(2,3) = z_min;
XYZ(3,1) = x_max;  XYZ(3,2) = y_min;  XYZ(3,3) = z_min;
XYZ(4,1) = x_max;  XYZ(4,2) = y_min;  XYZ(4,3) = z_max;
XYZ(5,1) = x_min;  XYZ(5,2) = y_min;  XYZ(5,3) = z_max;
XYZ(6,1) = x_min;  XYZ(6,2) = y_min;  XYZ(6,3) = z_min;
XYZ(7,1) = x_min;  XYZ(7,2) = y_max;  XYZ(7,3) = z_max;
XYZ(8,1) = x_max;  XYZ(8,2) = y_max;  XYZ(8,3) = z_max;

%refinement
x_max_ref = mesharg.ref_horiz /2;
x_min_ref = -x_max_ref;
y_max_ref =  x_max_ref;
y_min_ref =  x_min_ref;
z_max_ref =  0;
z_min_ref = -mesharg.ref_vert;
XYZ(9,1)  = x_min_ref; XYZ(9,2)  = y_max_ref; XYZ(9,3)  = z_min_ref;
XYZ(10,1) = x_max_ref; XYZ(10,2) = y_max_ref; XYZ(10,3) = z_min_ref;
XYZ(11,1) = x_max_ref; XYZ(11,2) = y_min_ref; XYZ(11,3) = z_min_ref;
XYZ(12,1) = x_max_ref; XYZ(12,2) = y_min_ref; XYZ(12,3) = z_max_ref;
XYZ(13,1) = x_min_ref; XYZ(13,2) = y_min_ref; XYZ(13,3) = z_max_ref;
XYZ(14,1) = x_min_ref; XYZ(14,2) = y_min_ref; XYZ(14,3) = z_min_ref;
XYZ(15,1) = x_min_ref; XYZ(15,2) = y_max_ref; XYZ(15,3) = z_max_ref;
XYZ(16,1) = x_max_ref; XYZ(16,2) = y_max_ref; XYZ(16,3) = z_max_ref;

% points list
pts_list = zeros(nb_point_box_edge,4);  %[nb x z y] NB: 'z' is in the screen plane in 2D and 'y' is depth (positive upward)
pts_list(:,1)  = (1:nb_point_box_edge)';
pts_list(:,2:end)  = XYZ; 

% facets list around the box domain
nb_face = 15;
face_4_box = zeros(nb_face,4+2);
face_4_box(1,:)  = [1  4  1 2 8 7 ]; % front
face_4_box(2,:)  = [2  4  2 3 4 8 ]; % rigth
face_4_box(3,:)  = [4  4  1 6 5 7 ]; % left
face_4_box(4,:)  = [6  4  1 2 3 6 ]; % bottom
face_4_box(5,:)  = [5  4  7 8 16 15 ]; % top front
face_4_box(6,:)  = [5  4  8 16 12 4 ]; % top right
face_4_box(7,:)  = [5  4  12 13 5 4 ]; % top back
face_4_box(8,:)  = [5  4  15 13 5 7 ]; % top left
face_4_box(9,:)  = [5  4  15 16 12 13 ]; % top center
face_4_box(10,:) = [11 4  9 10 16 15 ]; % front_intru
face_4_box(11,:) = [11 4  10 11 12 16 ]; % rigth intru
face_4_box(12,:) = [11 4  11 12 13 14 ]; % back intru
face_4_box(13,:) = [11 4  9 14 13 15 ]; % left intru
face_4_box(14,:) = [11 4  11 14 9 10 ]; % bottom intru
face_4_box(15,:) = [3  4  3 6 5 4 ]; % back
% 
% figure(1),clf,hold on
% for i = 1:nb_face 
%  plot3(pts_list(face_4_box(i,[3 4 5 6 3]),2),pts_list(face_4_box(i,[3 4 5 6 3]),3),pts_list(face_4_box(i,[3 4 5 6 3]),4),'-rd')
% end


%%%%%%%%%%%%%% GENERATE PTS AND FACE OF INTRUSION DOMAIN %%%%%%%%%%%%%%%%%%
b_el     = mesharg.ratio_ell; % Ellipse parameter
rad      = mesharg.lgth_ell/2 * b_el; % ellipse radius small axis

% Generate a base icosahedron mesh (Part of S2 Sampling Toolbox from Mathworks (Anton Semechko, 2015) )
ref_sub = mesharg.ref_sub_ellipse; %[1==42nodes , 2==162nodes, 3==642nodes, 4==2562nodes, 5==10242nodes, 6==40962nodes, 7==163842nodes, 8==655362nodes, 9==2621442nodes]
TR=IcosahedronMesh;
TR=SubdivideSphericalMesh(TR,ref_sub);
xyz = TR.X; 
tri = TR.Triangulation + nb_point_box_edge;
xyz(:,3) = (xyz(:,3) * mesharg.tck_ell/2) - mesharg.depth_ell ;  
xyz(:,[1 2]) = xyz(:,[1 2]) * rad;
figure(1),clf
h=trimesh(TR.Triangulation, xyz(:,1), xyz(:,2),xyz(:,3)); set(h,'EdgeColor','b','FaceColor','w')
axis equal
drawnow

% points list
nb_pts         = size(xyz,1);
pts_list_extra = [nb_point_box_edge+(1:nb_pts)' xyz];

% facet list
nb_face_3_ell = size(tri,1);
face_3_ell    = zeros(nb_face_3_ell ,3+2);
face_3_ell(:,1:2) = ones(nb_face_3_ell,1)*[7 3];
face_3_ell(:,3:end) = tri;


%%%%%%%%%%%%%%%%%%%%%   WRITING THE POLY FILE   %%%%%%%%%%%%%%%%%%%%%%%%%%%

% Points list
pts_list = [ pts_list ; pts_list_extra];
nb_pts   = size(pts_list,1) ;
% Face numbering
nb_face = size(face_4_box,1) + size(face_3_ell,1);
% nb_face = size(face_4_box,1) + size(face_3_sill,1) + size(face_3_topbot,1) ;
% Regions list
tmp = length(mesharg.max_tet_vol)+length(mesharg.ref_fact);
reg_list      = zeros(tmp, 6); % [nb x y z reg_nb reg_area]
reg_list(:,1) = (1:tmp)';
reg_list(1,2:end) = [0 0 -mesharg.depth_mod+1 1 mesharg.max_tet_vol];
reg_list(2,2:end) = [0 0 -1                   2 mesharg.max_tet_vol/mesharg.ref_fact];
% reg_list(3,2:end) = [0 0 -mesharg.depth_sill   2 mesharg.max_tet_vol];

% reg_list(end,:)   = [];
nb_reg = size(reg_list,1);
% Holes
hole_list =[1 0 0 -mesharg.depth_ell];
nb_hole   = 1;

% hole_list =[];
% nb_hole   = 0;

%%% WRITE INPUT TETGEN (.POLY FILE)
fid     = fopen([modelname,'.poly'],'w');
fprintf(fid,'%d 3 0 0\n', nb_pts);
fprintf(fid,'%d %15.12f %15.12f %15.12f\n', pts_list');
fprintf(fid,'%d 1\n', nb_face);
fprintf(fid,'1 0 %d\n%d %d %d %d %d\n', face_4_box');
tmp_str = ['1 0 %d\n', repmat('%d ',1, face_3_ell(1,2)+1 ),'\n'];
fprintf(fid, tmp_str, face_3_ell');
% tmp_str = ['1 0 %d\n', repmat('%d ',1, face_3_topbot(1,2)+1 ),'\n'];
% fprintf(fid, tmp_str, face_3_topbot');
fprintf(fid,'%d\n',nb_hole);
fprintf(fid,'%d %15.12f %15.12f %15.12f\n', hole_list');
fprintf(fid,'%d\n', nb_reg);
fprintf(fid,'%d %15.12f %15.12f %15.12f %d %d\n', reg_list');
fclose(fid);


%%% EXECUTE TETGEN
tetgen_path  = [SOLVER.path_libr,'/tetgen1.5.0/tetgen'];
tic
fprintf(1, '\nMeshing with Tetgen:  ');
[dum1, dum2] = system([tetgen_path,' -pqqAa ',modelname,'.poly']);
fprintf(1, [num2str(toc), ' sec.']);


%%% READ OUTPUT TETGEN
%NODES READING
fid =fopen(strcat(modelname,'.1.node'),'r');
tmp = fscanf(fid, '%d',4);
nnod = tmp(1);
MESH.NODES = fscanf(fid, '%e', [4, nnod]);
fclose(fid);
MESH.NODES(1,:)   = [];

%FACET READING
fid =fopen(strcat(modelname,'.1.face'),'r');
tmp = fscanf(fid, '%d',2);
nfac = tmp(1);
FACET = fscanf(fid, '%d', [5, nfac]);
fclose(fid);
MESH.node_markers = zeros(1,nnod);
indx_bfc = find(FACET(5,:)==1 | FACET(5,:)==2 | FACET(5,:)==16 | FACET(5,:)==32); % index boundary facet
indx_4   = find(FACET(5,:)==4); % bottom  model
indx_8   = find(FACET(5,:)==8); % surface model
MESH.node_markers(FACET(2:4,indx_bfc)) = repmat(FACET(5,indx_bfc),3,1);
MESH.node_markers(FACET(2:4,indx_4))   = repmat(FACET(5,indx_4),3,1);
MESH.node_markers(FACET(2:4,indx_8))   = repmat(FACET(5,indx_8),3,1);
MESH.FACES = FACET(2:end-1,:);
MESH.face_markers = FACET(end,:);

%ELEMS READING
fid        = fopen(strcat(modelname,'.1.ele'),'r');
tmp        = fscanf(fid, '%d',3);
nel        = tmp(1);
MESH.ELEMS = fscanf(fid, '%d',[6, nel]);
fclose(fid);
MESH.elem_markers = MESH.ELEMS(end,:);
MESH.ELEMS(1,:)   = [];
MESH.ELEMS(end,:) = [];
MESH.ELEMS	      = uint32(MESH.ELEMS);

%%% Deleting mesh files
delete([modelname,'.1.node'])
delete([modelname,'.1.face'])
delete([modelname,'.1.ele'])
delete([modelname,'.1.edge'])
delete([modelname,'.poly'])

%%% Mesh Visualization
%paraview_write('tet4', MESH.NODES,  MESH.ELEMS,  MESH.elem_markers,  MESH.node_markers,'toto', 0)


%%% Change Region numbering in case of inner/outer domain refinement 
MESH.elem_markers(MESH.elem_markers==2) = 1;

%%% FIND NEIGHBORGS IN MESH
if exist('MESH.NEIGHBORS','var')==0 
    tmp_elems = double(MESH.ELEMS(1:4,:));
    tmp_nodes = MESH.NODES(:,1:max(tmp_elems(:)));
    trep = triangulation(tmp_elems', tmp_nodes');
    MESH.NEIGHBORS = trep.neighbors()';
    MESH.NEIGHBORS(isnan(MESH.NEIGHBORS)) = 0;
    MESH.NEIGHBORS = uint32(MESH.NEIGHBORS);
end

%%% FIND NORMAL TO INTRUSION FACETS 
MESH = Normal_calc_intrusion(MESH, mesharg, SOLVER);

%%% update node markers based on facet numbering
indx_fb = find( MESH.face_markers==mesharg.face_id(1) | MESH.face_markers==mesharg.face_id(3) ); % index boundary facet (front and back)
indx_rl = find( MESH.face_markers==mesharg.face_id(2) | MESH.face_markers==mesharg.face_id(4) ); % index boundary facet (right and left)
indx_5  = find( MESH.face_markers==mesharg.face_id(5) ); % surface  model
indx_6  = find( MESH.face_markers==mesharg.face_id(6) ); % bottom model
indx_7  = find( MESH.face_markers==mesharg.face_id(7) ); % Intrusion body

MESH.node_markers = zeros(1,size(MESH.NODES,2));
nnpf = 3;  % for tet4 elements!
MESH.node_markers(MESH.FACES(:,indx_fb)) = repmat(MESH.face_markers(indx_fb),nnpf,1);
MESH.node_markers(MESH.FACES(:,indx_rl)) = repmat(MESH.face_markers(indx_rl),nnpf,1);
MESH.node_markers(MESH.FACES(:,indx_5))  = repmat(MESH.face_markers(indx_5),nnpf,1);
MESH.node_markers(MESH.FACES(:,indx_6))  = repmat(MESH.face_markers(indx_6),nnpf,1);
MESH.node_markers(MESH.FACES(:,indx_7))  = repmat(MESH.face_markers(indx_7),nnpf,1); 



%% REODERING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MRCM Reordering of NODES
spopts.symmetric = 0;
spopts.ndof      = 1;
K = sparse_create(MESH.ELEMS,1, spopts);
[perm,iperm] = mrcm(K);

MESH.NODES = MESH.NODES(:,perm);
MESH.node_markers = MESH.node_markers(:,perm);
MESH.FACES = reshape(iperm(MESH.FACES), size(MESH.FACES));
MESH.ELEMS = reshape(iperm(MESH.ELEMS), size(MESH.ELEMS));


% Geometric Reordering of ELEMS
[~, permel]  = sort(min(MESH.ELEMS)); 
% [~, permel]  = sort(mean(MESH_ep.ELEMS)); 
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

