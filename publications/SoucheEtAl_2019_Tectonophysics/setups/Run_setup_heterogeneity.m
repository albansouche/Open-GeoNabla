%% 2D ELASTO-PLASTIC MODEL TO STUDY HOST HETEROGENEOUS AROUND PRESSUIZED MAGMATIC INTRUSIONS:

% Part of Open-GeoNabla, copyright GPLv3, 2018
% https://github.com/albansouche/Open-GeoNabla/
% Physics of Geological Processes (PGP) , The NJORD Centre, Dept of Geosciences, University of Oslo
% Author: Alban Souche

% SCRIPT USED FOR THE PUBLICATION:
% Souche, A., Galland, O., Haug, Ø. T., & Dabrowski, M. (2019). Impact of 
% host rock heterogeneity on failure around pressurized conduits: Implications
% for finger-shaped magmatic intrusions. Tectonophysics.

% SHORT DETAILS ON SETUP AND SOLVER
% - Cylindrical / Finger shaped magmatic intrusion
% - Drucker-Prager Elasto-plasticity model 
% - Interface external force applied by increment to the host
% - Elasto-plasticity solved with Newton-Raphson algorithm with consistent tangent operator


% CLEARING ----------------------------------------------------------------
clear

% PATH TO LIBRARIES -------------------------------------------------------
restoredefaultpath
path_ext_libr = '../../../CODES/ext_libr';
path_2D_libr  = '../../../CODES/2D_libr';
addpath(genpath(path_ext_libr))
addpath(genpath(path_2D_libr))
warning off %make sure MUTILS is compiled with OpenMP for best performance

%--------------------------------------------------------------------------
% MESH ARGUMENTS ----------------------------------------------------------
%--------------------------------------------------------------------------

% geometry domain
MESHARG.X_mod = 100;  % [m]
MESHARG.Y_mod = 100;   % [m]
MESHARG.X_ref = 3.0;
MESHARG.Y_ref = 3.0;

% ellipstic intrusion (quater)
MESHARG.ell_tck    = 1; % [m]
MESHARG.ell_ratio  = 1;   % ratio between max and min ellispe axis
MESHARG.no_pts_tip = 200;  % nb of nods to discretize the crack tip

% refinement per area
MESHARG.max_tri_area = [ 1  0.00008 ]; % maximum triangle area [m2]

% segments id        [wall_b <%>; wall_r <%>; wall_t <%><->; wall_l <%>; wall_b_sill <.>; sill_contour <@>]
MESHARG.seg_id     = [1         ; 2         ; 3            ; 4         ; 5              ; 6               ]; % segments id (11 is forbidden)

% elements type solid meshes (better to keep 'tri7', alternative is 'tri3')
MESHARG.type_el_ep ='tri7';

%--------------------------------------------------------------------------
% PHYSICAL PARAMETERS -----------------------------------------------------
%--------------------------------------------------------------------------
PARAM.g         = 0*[0; -9.81]; % gravity of Earth [m/s2]
% elasto-plastic domain
PARAM.ep_rho    = 2700 ;  % density [kg/m3]
PARAM.ep_E      = 50e9 ;  % Young's modulus [Pa]
PARAM.ep_nu     = 0.25 ;  % Poisson's ratio
PARAM.ep_Co     =  1e7 ;  % Cohesion [Pa]
PARAM.ep_lambda = PARAM.ep_E.*PARAM.ep_nu./((1+PARAM.ep_nu).*(1-2*PARAM.ep_nu));
PARAM.ep_G      = PARAM.ep_E./(2*(1+PARAM.ep_nu));   % shear modulus
PARAM.ep_K      = PARAM.ep_E./(3*(1-2*PARAM.ep_nu)); % bulk modulus
PARAM.ep_H      = 0;           % isotropic hardening
PARAM.ep_phi    = 30*(pi/180); % frictional angle
PARAM.ep_xsi    = 30*(pi/180); % dilatancy angle (NB: param.xsi<=param.phi)
PARAM.ep_eta    = 3*tan(PARAM.ep_phi) / sqrt(9 + (12*tan(PARAM.ep_phi)*tan(PARAM.ep_phi) )); % DP/MC plane strain match
PARAM.ep_etabar = 3*tan(PARAM.ep_xsi) / sqrt(9 + (12*tan(PARAM.ep_xsi)*tan(PARAM.ep_xsi) )); % DP/MC plane strain match
PARAM.ep_xi     = 3 / sqrt(9 + (12*tan(PARAM.ep_phi)*tan(PARAM.ep_phi) ));                   % DP/MC plane strain match

%--------------------------------------------------------------------------
% BOUNDARY CONDITIONS -----------------------------------------------------
%--------------------------------------------------------------------------
BC.i_press   = 30e6; % [Pa] constant pressure imposed along boundary <.>
% elasto-plastic domain
BC.ep_wall_t = 'free_surf'; % 'free_surf' or 'free_slip' or 'no_slip'
BC.ep_wall_b = 'free_slip'; % 'free_surf' or 'free_slip' or 'no_slip'
BC.ep_wall_r = 'free_slip'; % 'free_surf' or 'free_slip' or 'no_slip'
BC.ep_wall_l = 'free_slip'; % 'free_surf' or 'free_slip' or 'no_slip'

%--------------------------------------------------------------------------
% SOLVER PARAMETERS -------------------------------------------------------
%--------------------------------------------------------------------------
SOLVER.ep_nb_incre    = 6;   % number of increment to apply external load to elasto-plastic domain
SOLVER.ep_rel_res_tol = 1e-9;  % relative error tolerance to elasto-plastic iterations
SOLVER.ep_NR          = 'NR_full_update'; % 'NR_full_update' or 'NR_line_search'
SOLVER.ep_max_NR_it   = 20;    % maximum number of NR iterations
SOLVER.nelblo         = 10000; % number of ele per block for matrix assemblage (MILAMIN approach)
SOLVER.nb_cpu         = 2;     % number of cpu to use for parallel calculations
SOLVER.solver_type    = 'direct_sym_chol'; % 'direct_sym_chol' / 'direct_sym' / 'direct_nonsym' / (develpment: 'pcg'/'AGMG')
if ~isequal(PARAM.ep_phi,PARAM.ep_xsi); SOLVER.solver_type = 'direct_nonsym'; end  % non-associativity == non-symmetric stiffness matrix

%--------------------------------------------------------------------------
% PLOTTING and SAVING (to ../resutls/)-------------------------------------
%--------------------------------------------------------------------------
set(0, 'DefaultFigureRenderer', 'zbuffer');
POST.folder_name = ['results/run_',num2str( max(0,length(dir('results/run*')))+1 )] ;
POST.fig_quality = '-r200';
POST.Hz_action   = 1; % frequency of executing the postprocessing
%             [ plot, save_plot, save_data ]     % 1 or 0
POST.action = [   0 ,     0    ,     0     ; ... % 'solid mesh'
                  0 ,     0    ,     0     ; ... % 'yield function'
                  0 ,     0    ,     0     ; ... % 'solid pressure'
                  0 ,     0    ,     0     ; ... % 'solid deviatoric stresses'
                  0 ,     0    ,     0     ; ... % 'solid total strain'
                  1 ,     0    ,     0     ; ... % 'solid 2nd invariant of dev. plastic/elastic strain'
                  0 ,     0    ,     0     ];    % 'solid displacement'

POST.xminp = 0; % xlimit of the plots
POST.xmaxp = 2.5; % xlimit of the plots
POST.yminp = 0; % ylimit of the plots
POST.ymaxp = 2.5; % ylimit of the plots


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% END OF USER SETUP -------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%% MESH GENERATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dimension
ndim          = 2;   % 2D simulation (plain strain formulation)
% generate meshes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MESH = generate_mesh_ell(MESHARG);

% mecanical mesh
MESH.EperN    = accumarray(MESH.ELEMS(:), ones(size(MESH.ELEMS(:))));
nnod          = size(MESH.NODES,2);
[nnodel,nel]  = size(MESH.ELEMS);
fprintf(1, ['\n Number of nodes:   ', num2str(nnod)]);
fprintf(1, ['\n Number of elems:   ', num2str(nel),'\n']);

%% Preprocessing
if  sum(sum(POST.action(:,2:3)))>0
    mkdir(POST.folder_name)
    save([POST.folder_name,'/data_input'])
end
% plotting fluid/mechanical meshes
if POST.action(1,1)==1
    figure(1),clf
    toplot = repmat(MESH.elem_markers,3,1);
    trisurf(reshape(1:3*nel,3, nel)', MESH.NODES(1,MESH.ELEMS(1:3,:)), MESH.NODES(2,MESH.ELEMS(1:3,:)), zeros(size(MESH.NODES(1,MESH.ELEMS(1:3,:)))), toplot);
    title('mechanical mesh'); view(2); axis image;
    drawnow
end

%% BOUNDARY CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bc on elasto-plastic domain
ind_seg_ux = [];
ind_seg_uy = [];
switch BC.ep_wall_t  % segments number 3
    case 'free_surf'
        % nothing to do
    case 'free_slip'
        % no ux constrain
        ind_seg_uy = [ ind_seg_uy  find(MESH.segment_markers==MESHARG.seg_id(3)) ];
    case 'no_slip'
        ind_seg_ux = [ ind_seg_ux  find(MESH.segment_markers==MESHARG.seg_id(3)) ];
        ind_seg_uy = [ ind_seg_uy  find(MESH.segment_markers==MESHARG.seg_id(3)) ];
end
switch BC.ep_wall_b  % segments number 1
    case 'free_surf'
        % nothing to do
    case 'free_slip'
        % no ux constrain
        ind_seg_uy = [ ind_seg_uy  find(MESH.segment_markers==MESHARG.seg_id(1)) ];
    case 'no_slip'
        ind_seg_ux = [ ind_seg_ux  find(MESH.segment_markers==MESHARG.seg_id(1)) ];
        ind_seg_uy = [ ind_seg_uy  find(MESH.segment_markers==MESHARG.seg_id(1)) ];
end
switch BC.ep_wall_r  % segments number 2
    case 'free_surf'
        % nothing to do
    case 'free_slip'
        % no uy constrain
        ind_seg_ux = [ ind_seg_ux  find(MESH.segment_markers==MESHARG.seg_id(2)) ];
    case 'no_slip'
        ind_seg_ux = [ ind_seg_ux  find(MESH.segment_markers==MESHARG.seg_id(2)) ];
        ind_seg_uy = [ ind_seg_uy  find(MESH.segment_markers==MESHARG.seg_id(2)) ];
end
switch BC.ep_wall_l  % segments number 4
    case 'free_surf'
        % nothing to do
    case 'free_slip'
        % no uy constrain
        ind_seg_ux = [ ind_seg_ux  find(MESH.segment_markers==MESHARG.seg_id(4)) ];
    case 'no_slip'
        ind_seg_ux = [ ind_seg_ux  find(MESH.segment_markers==MESHARG.seg_id(4)) ];
        ind_seg_uy = [ ind_seg_uy  find(MESH.segment_markers==MESHARG.seg_id(4)) ];
end
BC.ep_bc_nod_ux  = unique(MESH.SEGMENTS(:,ind_seg_ux));
BC.ep_bc_val_ux  = zeros(size(BC.ep_bc_nod_ux));
BC.ep_bc_nod_uy  = unique(MESH.SEGMENTS(:,ind_seg_uy));
BC.ep_bc_val_uy  = zeros(size(BC.ep_bc_nod_uy));
BC.ep_ind = [ndim*(BC.ep_bc_nod_ux(:)-1)+1 ; ndim*(BC.ep_bc_nod_uy(:)-1)+2];
BC.ep_val = [BC.ep_bc_val_ux(:) ; BC.ep_bc_val_uy(:)];

%% Integration and evaluation points
if nnodel==3
    SOLVER.nip         = 1;
    SOLVER.neva        = 1; % Constant discontinuous strain field
    SOLVER.pts_eva_loc = [0,0; 1,0; 0,1];
elseif nnodel==7
    SOLVER.nip         = 6;
    SOLVER.neva        = 6; % Linear discontinuous strain field
    SOLVER.pts_eva_loc = [0,0; 1,0; 0,1; .5,.5; 0,.5; .5,0; 1/3,1/3];
end
[ieuv, ~] = ip_triangle(SOLVER.neva);


% Traction BC (total external force vector 'Fext')
MESH.normal(:,MESH.segment_markers<MESHARG.seg_id(6))=0; % normal vectors only at intrusion
SIGMA = zeros(3,nnod);
tmp_nods = unique(MESH.SEGMENTS(:,MESH.segment_markers==6));
SIGMA(1,tmp_nods) = -BC.i_press;
SIGMA(2,tmp_nods) = -BC.i_press;
seg_id        = MESHARG.seg_id(6);
seg_id_corner = [MESHARG.seg_id(1) MESHARG.seg_id(4)];
Fext   = ContourTraction(MESH, SIGMA, seg_id, seg_id_corner); % Rhs Traction vector (external forces)
figure(1),hold on,quiver(MESH.NODES(1,:),MESH.NODES(2,:),Fext(1:2:end)',Fext(2:2:end)','r' )

%% Solve lithostatic pressure
P_litho = initial_litho_zlevel(MESH, PARAM, SOLVER, POST);


%% Adding pertubation
load('isotropic_perturbation_lambda1.00.mat') % loading the perturbation map
noise_interp_mesh =  interp2( X, Y, Zf, MESH.NODES(1,MESH.ELEMS(7,:)), MESH.NODES(2,MESH.ELEMS(7,:)) );
tmp = noise_interp_mesh;
tmp(isnan(noise_interp_mesh)) = [];
std_n = std(tmp)
noise_interp_mesh = (1/std_n)*noise_interp_mesh; %rescaling by std of sampled mesh values

%figure(22),clf
%toplot = repmat(noise_interp_mesh,3,1);
%trisurf(reshape(1:3*size(MESH.ELEMS,2),3, size(MESH.ELEMS,2))', MESH.NODES(1,MESH.ELEMS(1:3,:)), MESH.NODES(2,MESH.ELEMS(1:3,:)), zeros(size(MESH.NODES(1,MESH.ELEMS(1:3,:)))), toplot);
%view(2)

PARAM.noise_amp    = 0.05; % relative to the cohesion (0.05 == 10% signal variation of the standard deviation)
PARAM.noise = PARAM.ep_Co(MESH.elem_markers) .* PARAM.noise_amp .* noise_interp_mesh;
PARAM.noise(isnan(PARAM.noise)) = 0;

%% Solve the EP problem
U_ep    = EP_DruckerPrager_calc_dSdU_ell(MESH, PARAM, BC, P_litho, Fext, SOLVER, POST);


