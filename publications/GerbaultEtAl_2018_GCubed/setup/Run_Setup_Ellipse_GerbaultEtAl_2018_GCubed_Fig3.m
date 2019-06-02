
%% SCRIPT TO EXECUTE 3D CALCULATIONS USED IN THE PUBLICATION:
% Gerbault et al., 2018, Three dimensional failure patterns around an inflating magmatic chamber. Geochemistry, Geophysics, Geosystems  

% Part of Open-GeoNabla, copyright GPLv3, 2018
% https://github.com/albansouche/Open-GeoNabla/
% Physics of Geological Processes (PGP) , The NJORD Centre, Dept of Geosciences, University of Oslo
% Author: Alban Souche

%%  OPENGEONABLA ----------------------------------------------------------
%   Elasto-Plasticity (Drucker-Prager)
%   3D finite element formulation (fully implicit)
%   Code build partially upon MILAMIN 1.0.1, Mutils 0.4.2, Tetgen1.5.0, FOLDER 1.0
%   Elasto-plastic implementation inspired from HYPLAS_v2.0
% -------------------------------------------------------------------------
%   DETAILS:
%   - Drucker-Prager failure criteria
%   - Incremental loading procedure to the external forces
%   - Elasto-plasticity solved with Newton-Raphson algorithm with consistent tangent operator
% -------------------------------------------------------------------------


% CLEARING ----------------------------------------------------------------
clear

% PATH TO LIBRARIES -------------------------------------------------------
path_ext_libr = fullfile('..','..','..','CODES','ext_libr');
path_3D_libr  = fullfile('..','..','..','CODES','3D_libr');
addpath(genpath(path_ext_libr))
addpath(genpath(path_3D_libr))

%--------------------------------------------------------------------------
% MESH ARGUMENTS ----------------------------------------------------------
%--------------------------------------------------------------------------
iR    = 500;  % Radius ellipse [m]
iD    = 2000; % Depth ellipse [m]
iLamb = linspace(0.2,3,20); % Ratio (lambda) Prolate/Oblate

iRes  = [5; 1e2]; % for  MESHARG.no_pts_ell and MESHARG.ref_fact
                  % iRes(1)=5 == LowRes ; iRes(1)=7 == HiRes

Data_Max_Uplift = zeros(length(iR),length(iD),length(iLamb),size(iRes,2));
Data_dPc        = zeros(length(iR),length(iD),length(iLamb),size(iRes,2));
Data_dPc_scaled = zeros(length(iR),length(iD),length(iLamb),size(iRes,2));
Data_dPc_ini    = zeros(length(iR),length(iD),length(iLamb),size(iRes,2));

% Example of Fig3b (Gerbault et al., GCube, 2018)
figure(11),clf
% DATA from Muriel  %%%%%%%%%%%%%%%
cw_la=[0.25 0.5 0.6  0.80 1.0 1.20 1.45 1.65 1.85 2 2.1 2.4 2.5 2.75 3];
cw_DP=[ 22  26  32   38   44   48  44  38  34  30 26  24  20.5  18 18];
plot ( cw_la, cw_DP,'b-*'); hold on;
mu1_la=[ 0.2    0.6    1  1.5    2  3.];
mu1_DP=[ 32.5  32.5   45  42    35  20];
plot (mu1_la,mu1_DP,'g-*'); hold on;
xlabel('lambda shape ratio');
ylabel('DP wall failure');
legend('CW r=0.5km','Adeli r=0.5km')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for idepth = 1:length(iD)
    
    for irad = 1:length(iR)
        
        for ilamb = 1:length(iLamb)
            
            for ires = 1:size(iRes,2)
                
                % geometry domain and layer
                MESHARG.width_mod  = 160e3; % [m]
                MESHARG.depth_mod  = 80e3; % [m]
                MESHARG.depth_ell = iD(idepth);  % [m]
                MESHARG.tck_ell   = 2*iR(irad);  % [m]
                MESHARG.lgth_ell  = 2*iR(irad)*iLamb(ilamb);  % [m]
                
                % sill intrusion
                MESHARG.ratio_ell  = 1;  % ellipse param.
                MESHARG.ref_sub_ellipse = iRes(1,ires); % refinement along ellipse interface
                
                % facets id       [wall_f; wall_r; wall_b; wall_l; wall_top; wall_bot; intrusion_contour]
                MESHARG.face_id = [1     ; 2     ; 3     ; 4     ; 5       ; 6       ; 7                ]; % facet id (11 is forbidden)
                
                % element refinement per area
                MESHARG.max_tet_vol = 2e10; % maximum triangle area [m2]
                
                % sepcific refinement inner/outer domain (comment the 2 lines if not needed)
                MESHARG.ref_horiz = 20000; %2 * MESHARG.lgth_sill * MESHARG.sill_y_ell;
                MESHARG.ref_vert  = 10000; %2 * MESHARG.depth_sill;
                MESHARG.ref_fact  = iRes(2,ires); % ref factor between inner/outer domains
                
                
                % elements type of fluid and solid meshes (better to keep 'tri7', alternative is 'tri3')
                MESHARG.type_el_ep ='tet4'; % Only 'tet4'
                
                %--------------------------------------------------------------------------
                % PHYSICAL PARAMETERS -----------------------------------------------------
                %--------------------------------------------------------------------------
                PARAM.g         = [0; 0 ; -9.81]; % gravity of Earth [m/s2]
                PARAM.ep_rho    = 2600 ;  % density [kg/m3]
                PARAM.ep_E      = 75e9 ;  % Young's modulus [Pa]
                PARAM.ep_nu     = 0.25 ;  % Poisson's ratio
                PARAM.ep_Co     =  5e6 ;  % Cohesion [Pa]
                PARAM.ep_lambda = PARAM.ep_E.*PARAM.ep_nu./((1+PARAM.ep_nu).*(1-2*PARAM.ep_nu)); % (1st Lam? coeff)
                PARAM.ep_G      = PARAM.ep_E./(2*(1+PARAM.ep_nu));   % shear modulus (==2nd Lam? coeff)
                PARAM.ep_K      = PARAM.ep_E./(3*(1-2*PARAM.ep_nu)); % bulk modulus
                PARAM.ep_H      = 0;           % isotropic hardening
                PARAM.ep_phi    = 35*(pi/180); % frictional angle
                PARAM.ep_xsi    = 35*(pi/180); % dilatancy angle (NB: associative plasticity xsi=phi)
                % Changed to match CurrentiWilliams_2014
                PARAM.ep_eta    = 6*sin(PARAM.ep_phi) / ( sqrt(3) * (3-sin(PARAM.ep_phi) ));
                PARAM.ep_etabar = 6*sin(PARAM.ep_xsi) / ( sqrt(3) * (3-sin(PARAM.ep_xsi) ));
                PARAM.ep_xi     = 6*cos(PARAM.ep_phi) / ( sqrt(3) * (3-sin(PARAM.ep_phi) ));
                
                %--------------------------------------------------------------------------
                % BOUNDARY CONDITIONS -----------------------------------------------------
                %--------------------------------------------------------------------------
                BC.i_press    = 300e6; % [Pa] constant pressure imposed along boundary <.>
                % elasto-plastic domain
                BC.ep_wall_t  = 'free_surf'; % top, 'free_surf' or 'free_slip' or 'no_slip'
                BC.ep_wall_b  = 'free_slip'; % bottom, 'free_surf' or 'free_slip' or 'no_slip'
                BC.ep_wall_rl = 'free_slip'; % lateral right and left walls , 'free_surf' or 'free_slip' or 'no_slip'
                BC.ep_wall_fb = 'free_slip'; % lateral front and back walls , 'free_surf' or 'free_slip' or 'no_slip'
                
                %--------------------------------------------------------------------------
                % SOLVER PARAMETERS -------------------------------------------------------
                %--------------------------------------------------------------------------
                SOLVER.ep_nb_incre    = BC.i_press/1e6;   % number of increment to apply external load to elasto-plastic domain
                SOLVER.ep_rel_res_tol = 1e-11;  % relative error tolerance to elasto-plastic iterations
                SOLVER.ep_NR          = 'NR_full_update'; % 'NewtonRaphson full update'
                SOLVER.ep_max_NR_it   = 100;    % maximum number of NR iterations
                SOLVER.nelblo         = 10000; % number of elem per block for matrix assemblage (MILAMIN approach)
                SOLVER.nb_cpu         = 1;     % number of cpu to use for parallel calculations
                SOLVER.solver_type    = 'PCG'; 
                
                %--------------------------------------------------------------------------
                % PLOTTING and SAVING (to ../resutls/)-------------------------------------
                %--------------------------------------------------------------------------
                set(0, 'DefaultFigureRenderer', 'zbuffer');
                POST.folder_name = ['../results/Ell_Intru_',num2str( max(1,length(dir('../results/Ell_Intru*'))+1) )] ;
                POST.Hz_action   = 10; % frequency of executing the postprocessing
                %             [ save_data ]  % 1 or 0
                POST.action = [ 0 ; ... % 'mesh'
                                0 ; ... % 'yield function'
                                0 ; ... % 'pressure'
                                0 ; ... % 'deviatoric stresses (second invariant)'
                                0 ; ... % 'deviatoric strain (second invariant)'
                                0 ; ... % 'volumetric strain'
                                0 ];    % 'displacement'
                
                %--------------------------------------------------------------------------
                %--------------------------------------------------------------------------
                % END OF USER SETUP -------------------------------------------------------
                %--------------------------------------------------------------------------
                %--------------------------------------------------------------------------
                
                %% MESH GENERATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % dimension
                ndim          = 3;   % 3D simulation
                % generate meshes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                modelname = [pwd,'/mesh_ell'];
                SOLVER.path_libr = path_ext_libr; 
                MESH = generate3Dmesh_ellipsoid(MESHARG, SOLVER, modelname);
                MESH.EperN    = accumarray(MESH.ELEMS(:), ones(size(MESH.ELEMS(:))));

                nnod          = size(MESH.NODES,2);
                [nnodel,nel]  = size(MESH.ELEMS);
                fprintf(1, ['\nNumber of nodes    :  ', num2str(nnod)]);
                fprintf(1, ['\nNumber of elems    :  ', num2str(nel),'\n\n']);
                
                            %% Preprocessing
                            if  sum(POST.action(:))>0
                                mkdir(POST.folder_name)
                                save([POST.folder_name,'/data_input'])
                            end
                            % plotting fluid/mechanical meshes
                            if POST.action(1)==1
                                paraview_write('tet4', MESH.NODES,  MESH.ELEMS,  MESH.elem_markers,  MESH.node_markers, [POST.folder_name,'/mesh_ini'], 0)
                            end
                
                            % bc on elasto-plastic domain
                            ind_node_ux = [];
                            ind_node_uy = [];
                            ind_node_uz = [];
                            switch BC.ep_wall_t  % nodes on facet number 5
                                case 'free_surf' % nothing to do
                                case 'free_slip' % no ux and uy constrain
                                    ind_tmp     = find(MESH.node_markers==5);
                                    ind_node_uz = [ ind_node_uz  ind_tmp ];
                                case 'no_slip'   % constrain in all direction
                                    ind_tmp     = find(MESH.node_markers==5);
                                    ind_node_ux = [ ind_node_ux  ind_tmp ];
                                    ind_node_uy = [ ind_node_uy  ind_tmp ];
                                    ind_node_uz = [ ind_node_uz  ind_tmp ];
                            end
                            switch BC.ep_wall_b  % nodes on facet number 6
                                case 'free_surf' % nothing to do
                                case 'free_slip' % no ux uy constrain
                                    ind_tmp     = find(MESH.node_markers==6);
                                    ind_node_uz = [ ind_node_uz  ind_tmp ];
                                case 'no_slip'   % constrain in all direction
                                    ind_tmp     = find(MESH.node_markers==6);
                                    ind_node_ux = [ ind_node_ux  ind_tmp ];
                                    ind_node_uy = [ ind_node_uy  ind_tmp ];
                                    ind_node_uz = [ ind_node_uz  ind_tmp ];
                            end
                            switch BC.ep_wall_fb  % nodes on facet number 1 and 3 (front and back)
                                case 'free_surf'  % nothing to do
                                case 'free_slip'  % no uy and uz constrain
                                    ind_tmp     = find(or(MESH.node_markers==1, MESH.node_markers==3));
                                    ind_node_ux = [ ind_node_ux  ind_tmp ];
                                case 'no_slip'    % constrain in all direction
                                    ind_tmp     = find(or(MESH.node_markers==1, MESH.node_markers==3));
                                    ind_node_ux = [ ind_node_ux  ind_tmp ];
                                    ind_node_uy = [ ind_node_uy  ind_tmp ];
                                    ind_node_uz = [ ind_node_uz  ind_tmp ];
                            end
                            switch BC.ep_wall_rl  % nodes on facet number 2 and 4 (right and left)
                                case 'free_surf'  % nothing to do
                                case 'free_slip'  % no ux and uz constrain
                                    ind_tmp     = find(or(MESH.node_markers==2, MESH.node_markers==4));
                                    ind_node_uy = [ ind_node_uy  ind_tmp ];
                                case 'no_slip'    % constrain in all direction
                                    ind_tmp     = find(or(MESH.node_markers==2, MESH.node_markers==4));
                                    ind_node_ux = [ ind_node_ux  ind_tmp ];
                                    ind_node_uy = [ ind_node_uy  ind_tmp ];
                                    ind_node_uz = [ ind_node_uz  ind_tmp ];
                            end
                            BC.ep_bc_nod_ux  = unique(ind_node_ux);
                            BC.ep_bc_val_ux  = zeros(size(BC.ep_bc_nod_ux));
                            BC.ep_bc_nod_uy  = unique(ind_node_uy);
                            BC.ep_bc_val_uy  = zeros(size(BC.ep_bc_nod_uy));
                            BC.ep_bc_nod_uz  = unique(ind_node_uz);
                            BC.ep_bc_val_uz  = zeros(size(BC.ep_bc_nod_uz));
                            BC.ep_ind = [ ndim*(BC.ep_bc_nod_ux(:)-1)+1 ; ndim*(BC.ep_bc_nod_uy(:)-1)+2 ; ndim*(BC.ep_bc_nod_uz(:)-1)+3 ];
                            BC.ep_val = [ BC.ep_bc_val_ux(:)            ; BC.ep_bc_val_uy(:)            ; BC.ep_bc_val_uz(:)            ];
                
                
                            %% Integration and evaluation points
                            SOLVER.nip         = 1;
                            SOLVER.neva        = 1; % Constant discontinuous strain field
                            SOLVER.pts_eva_loc = [0,0,0; 1,0,0; 0,1,0; 0,0,1];
                            [ieuv, ~]          = ip_tetra(SOLVER.neva);
                
                            %% Traction BC (total external force vector 'Fext')
                            SIGMA = zeros(6,nnod);
                            i_nodes = find(MESH.node_markers==7);
                            SIGMA(1,i_nodes) = - BC.i_press; % xx
                            SIGMA(2,i_nodes) = - BC.i_press; % yy
                            SIGMA(3,i_nodes) = - BC.i_press; % zz
                            face_id = 7; % traction face_id 
                            Fext  = ContourTraction(MESH, SIGMA, face_id); % Rhs Traction vector (external forces)

                            %% Mechanical solver
                            P_litho = initial_litho(MESH, PARAM, SOLVER, POST);
                            [U_ep, dP_c_i, dP_c] = EP_DruckerPrager_Critic_dP(MESH, PARAM, BC, P_litho, Fext, SOLVER, POST);
                
                            Data_dPc_ini(irad,idepth,ilamb,ires)    = dP_c_i*1e6;
                            Data_dPc(irad,idepth,ilamb,ires)        = dP_c*1e6;
                            Data_dPc_scaled(irad,idepth,ilamb,ires) = dP_c*1e6/abs(PARAM.g(3)*PARAM.ep_rho*iD(idepth));
                            
                            uz = U_ep(3:3:end);
                            indx = unique(find(MESH.node_markers==5));
                            Data_Max_Uplift(irad,idepth,ilamb,ires) = max(uz(indx));
                            
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            figure(11),hold on
                            plot(iLamb(ilamb),dP_c, '-ko','MarkerFaceColor','r','MarkerSize',8)
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
            end
            
        end
        
    end
end


% save('Data_Fig3','Data_Max_Uplift','Data_dPc_scaled','Data_dPc','Data_dPc_ini','iR','iD','iLamb','iRes')


figure(11),clf, hold on

% DATA from Muriel  %%%%%%%%%%%%%%%
cw_la=[0.25 0.5 0.6  0.80 1.0 1.20 1.45 1.65 1.85 2 2.1 2.4 2.5 2.75 3];
cw_DP=[ 22  26  32   38   44   48  44  38  34  30 26  24  20.5  18 18];
plot ( cw_la, cw_DP,'b-*'); hold on;

mu1_la=[ 0.2    0.6    1  1.5    2  3.];
mu1_DP=[ 32.5  32.5   45  42    35  20];
plot (mu1_la,mu1_DP,'g-*'); hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(iLamb,squeeze(Data_dPc(1,1,:,1))/1e6, '-ko','MarkerFaceColor','r','MarkerSize',8)
xlabel('lambda shape ratio')
ylabel('DP wall failure')
legend('CW r=0.5km','Adeli r=0.5km','Adhoc r=0.5km (lower res)','Adhoc r=0.5km (higer res)')






