%% SCRIPT TO EXECUTE 2D ELASTIC BENCHMARK OF A PRESSUIZED MAGMA RESERVOIR :

% Part of Open-GeoNabla, copyright GPLv3, 2018
% https://github.com/albansouche/Open-GeoNabla/
% Physics of Geological Processes (PGP) , The NJORD Centre, Dept of Geosciences, University of Oslo
% Author: Alban Souche

% SHORT DETAILS ON SETUP AND SOLVER
% - Cylindrical magma reservoir
% - Elastic solver

% Solution from:
% Gerbault, M., F. Cappa, and R. Hassani (2012), Correction to "Elasto-plastic and hydromechanical models of failure around an
% infinitely long magma chamber", Geochem. Geophys. Geosyst., 13, Q10015, doi:10.1029/2012GC004464.


%% CLEARING ----------------------------------------------------------------
clear

% PATH TO LIBRARIES -------------------------------------------------------
restoredefaultpath
path_ext_libr = '../../../CODES/ext_libr';
path_2D_libr  = '../../../CODES/2D_libr';
addpath(genpath(path_ext_libr))
addpath(genpath(path_2D_libr))

% Coefficient recaling the length and width of the domain
Lres = logspace(8,1,20);
Lres = fliplr(log10(Lres));

% Coefficient recaling the space discretisation of the domain
Res = .5;

erx =  zeros(1,length(Lres));
ery =  zeros(1,length(Lres));

% Loop over different lengths and widths of the domain
for ii = 1:length(Lres)
    
    % Loop over different space discretisation resolution
    for i = 1:length(Res)
        
        %--------------------------------------------------------------------------
        % MESH ARGUMENTS ----------------------------------------------------------
        %--------------------------------------------------------------------------
        
        %     %''''''''' Free surface ''''''''%
        %     %                               %   <'>  free surface
        %     %-------------                  %   <%>  displacement boundary conditions for elasto-plastic domain
        %       .           |                 %   <.>  boundary where external load is applied
        %          .        |                 %   <(number)> region for mesh refinement (reset to 1 everywhere afterwards)
        %            .      |                 %
        %      over   .     |                 %
        %   pressured .     |                 %
        %   reservoir .     |                 %
        %            .      |                 %
        %          .        |                 %
        %       .       (2) |  (1)            %
        %     %-------------                  %
        %     %                               %
        %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % geometry domain
        MESHARG.Lx           = Lres(ii)*100e3; % depth model [m]
        MESHARG.Ly           = Lres(ii)*100e3;  % width model
        MESHARG.Mag_depth    = 7e3;
        MESHARG.Mag_rad      = 2e3;  % radius of the crack tip
        MESHARG.Mag_nb_pts   = 1000;  % nb of nods along the crack tip
        MESHARG.x_ref        = 8e3 + MESHARG.Mag_rad;   % refined depth model
        MESHARG.y_ref        = 2e3 + MESHARG.Mag_depth + MESHARG.Mag_rad; % refined width model
        % segments id          [wall_b; wall_r; wall_t; wall_l; wall_l_sill; sill_contour]
        MESHARG.seg_id       = [     1;      2;      3;      4;           5;            6]; % segments id
        MESHARG.seg_dx_max   = [    50;     50;     50;     50;          10;          100]*1e4; % max dx along segments
        MESHARG.max_tri_area = Res(i)*[1e7, 1e5]; % maximum triangle area [area1, area2]
        MESHARG.type_el_ep   ='tri7';
        
        %--------------------------------------------------------------------------
        % PHYSICAL PARAMETERS -----------------------------------------------------
        %--------------------------------------------------------------------------
        PARAM.g         = [0; -9.81]; % gravity of Earth [m/s2]
        PARAM.ep_rho    = 2500;      % density [kg/m3]
        PARAM.ep_E      = 50*1e9;    % Young's modulus [Pa]
        PARAM.ep_nu     = 0.25 ;     % Poisson's ratio
        PARAM.ep_lambda = PARAM.ep_E.*PARAM.ep_nu./((1+PARAM.ep_nu).*(1-2*PARAM.ep_nu));
        PARAM.ep_G      = PARAM.ep_E./(2*(1+PARAM.ep_nu));   % shear modulus
        PARAM.ep_K      = PARAM.ep_E./(3*(1-2*PARAM.ep_nu)); % bulk modulus
        PARAM.ep_phi    = nan;
        PARAM.ep_xsi    = nan;
        %--------------------------------------------------------------------------
        % BOUNDARY CONDITIONS -----------------------------------------------------
        %--------------------------------------------------------------------------
        % viscous domain  (pressure bc)
        BC.i_press   = 50e6; % [Pa] constant pressure imposed along boundary <.>
        % elasto-plastic domain
        BC.ep_wall_t = 'free_surf'; % 'free_surf' or 'free_slip' or 'no_slip'
        BC.ep_wall_b = 'free_slip'; % 'free_surf' or 'free_slip' or 'no_slip'
        BC.ep_wall_r = 'no_slip'; % 'free_surf' or 'free_slip' or 'no_slip'
        BC.ep_wall_l = 'free_slip'; % 'free_surf' or 'free_slip' or 'no_slip'
        
        %--------------------------------------------------------------------------
        % SOLVER PARAMETERS -------------------------------------------------------
        %--------------------------------------------------------------------------
        SOLVER.nelblo      = 10000; % number of ele per block for matrix assemblage (MILAMIN approach)
        SOLVER.nb_cpu      = 2;     % number of cpu to use for parallel calculations
        SOLVER.solver_type = 'direct_sym_chol'; % 'direct_sym_chol' / 'direct_sym' / 'direct_nonsym' / (develpment: 'pcg'/'AGMG')
        %--------------------------------------------------------------------------
        % PLOTTING and SAVING (to ../resutls/)-------------------------------------
        %--------------------------------------------------------------------------
        set(0, 'DefaultFigureRenderer', 'zbuffer');
        POST.folder_name = ['../results/Gerbault2012_',num2str( max(1,length(dir('../results/Gerbault2012*'))+1) )] ;
        POST.fig_quality = '-r200';
        POST.Hz_action   = 1; % frequency of executing the postprocessing
        %             [ plot, save_plot, save_data ]     % 1 or 0
        POST.action = [   0 ,     0    ,     0     ; ... % 'fluid solid meshes'
                          0 ,     0    ,     0     ; ... % 'fluid velocity'
                          0 ,     0    ,     0     ; ... % 'fluid pressure'
                          0 ,     0    ,     0     ; ... % 'yield function'
                          0 ,     0    ,     0     ; ... % 'solid pressure'
                          0 ,     0    ,     0     ; ... % 'solid deviatoric stresses'
                          0 ,     0    ,     0     ; ... % 'solid deviatoric strain'
                          0 ,     0    ,     0     ; ... % 'solid volumetric strain'
                          0 ,     0    ,     0     ; ... % 'solid accumulated plastic strain'
                          0 ,     0    ,     0     ];    % 'solid displacement'
        POST.xminp = 0; % xlimit of the plots
        POST.xmaxp = 14e3; % xlimit of the plots
        POST.yminp = -14e3; % ylimit of the plots
        POST.ymaxp = 0; % ylimit of the plots
      
        
        %--------------------------------------------------------------------------
        %--------------------------------------------------------------------------
        % END OF USER SETUP -------------------------------------------------------
        %--------------------------------------------------------------------------
        %--------------------------------------------------------------------------
        
        
        %% MESH GENERATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % dimension
        ndim          = 2;   % 2D simulation (plain strain formulation)
        % generate mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        MESH = generate_mesh_Gerbault2012(MESHARG);
        % mecanical mesh reordering and more
        [MESH,~,~]    = mesh_reorder_amd(MESH);
        MESH.EperN    = accumarray(MESH.ELEMS(:), ones(size(MESH.ELEMS(:))));
        nnod          = size(MESH.NODES,2);
        [nnodel,nel]  = size(MESH.ELEMS);
        MESH  = seg_manip(MESH, MESHARG.type_el_ep, SOLVER);
        MESH  = NormalEdges(MESH);
        fprintf(1, ['\n Number of nodes:   ', num2str(nnod)]);
        fprintf(1, ['\n Number of elems:   ', num2str(nel),'\n']);
        
        
        %% BOUNDARY CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        
        % bc from over pressured reservoir
        BC.res_inter = unique(MESH.SEGMENTS(:, MESH.segment_markers==MESHARG.seg_id(6))); % dyke inlet (segments 5)
        SIGMA = zeros(3,nnod);
        SIGMA(1,BC.res_inter) = BC.i_press;
        SIGMA(2,BC.res_inter) = BC.i_press;
        seg_id        = MESHARG.seg_id(6);
        seg_id_corner = MESHARG.seg_id(4);
        Fext  = - ContourTraction(MESH, SIGMA, seg_id, seg_id_corner); % Rhs Traction vector (external forces)
        
        
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
        
        %% Mechanical solver
        P_litho = initial_litho(MESH, PARAM, SOLVER, POST);
        
        %% Boundary conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        BC_e.ind = BC.ep_ind;   % for elastic predictor increment
        BC_e.val = BC.ep_val;
        
        %% Solve elastic strain U %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Elastic tangent matix
        indx_p  = false(SOLVER.nip,nel);
        D_el_ip = Tangent_Elastoplastic_DruckerPrager(MESH, indx_p, [], [], [], PARAM, SOLVER.nip);
        % Assemble stiffness matrix and rhs
        elast = 1; % '1'== elastic only ; '0'==elastoplastic
        [A_all, ~, InvJ_all, DetJ_all] = matrix_entries_EP(MESH, D_el_ip, PARAM, SOLVER, elast);
        %Solve  elastic problem
        U = mechanical_solver(MESH, A_all, BC_e, Fext, SOLVER); % Fext Traction BC
        
        
        %% POST PROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        xminp = POST.xminp;
        xmaxp = POST.xmaxp;
        yminp = POST.yminp;
        ymaxp = POST.ymaxp;
        
        
        figure(1),clf
        subplot(2,2,1)
        toplot = U(1:2:end-1);
        toplot = toplot(MESH.ELEMS(1:3,:));
        trisurf(reshape(1:3*nel,3, nel)', MESH.NODES(1,MESH.ELEMS(1:3,:)), MESH.NODES(2,MESH.ELEMS(1:3,:)), zeros(size(MESH.NODES(1,MESH.ELEMS(1:3,:)))), toplot);
        title('Ux'); view(2); colorbar; shading interp; grid on; box on; axis equal
        xlabel('Distance [m]','FontSize',10,'FontWeight','bold')
        ylabel('Depth [m]' ,'FontSize',10,'FontWeight','bold')
        subplot(2,2,2)
        toplot = U(2:2:end);
        toplot = toplot(MESH.ELEMS(1:3,:));
        trisurf(reshape(1:3*nel,3, nel)', MESH.NODES(1,MESH.ELEMS(1:3,:)), MESH.NODES(2,MESH.ELEMS(1:3,:)), zeros(size(MESH.NODES(1,MESH.ELEMS(1:3,:)))), toplot);
        title('Uy'); view(2); colorbar; shading interp; grid on; box on; axis equal
        xlabel('Distance [m]','FontSize',10,'FontWeight','bold')
        ylabel('Depth [m]' ,'FontSize',10,'FontWeight','bold')
        subplot(2,2,3), hold on
        surf_nod = unique(MESH.SEGMENTS(:, find(MESH.segment_markers==3)));
        UX = U(1:2:end);
        plot(MESH.NODES(1,surf_nod),UX(surf_nod),'kd')
        plot([min(MESH.NODES(1,:)) max(MESH.NODES(1,:))] , [0 0],'-k')
        %xlim([xminp xmaxp])
        grid on, box on
        title('Surface horizontal displacement [m]')
        xlabel('Distance [m]','FontSize',10,'FontWeight','bold')
        subplot(2,2,4), hold on
        surf_nod = unique(MESH.SEGMENTS(:, find(MESH.segment_markers==3)));
        UY = U(2:2:end);
        plot(MESH.NODES(1,surf_nod),UY(surf_nod),'kd')
        plot([min(MESH.NODES(1,:)) max(MESH.NODES(1,:))] , [0 0],'-k')
        %xlim([xminp xmaxp])
        grid on, box on
        title('Surface vertical displacement [m]')
        xlabel('Distance [m]','FontSize',10,'FontWeight','bold')
        drawnow
        
        
        %%%%%% PLOT NUMERICAL SOLUTION %%%%%%%%%%%%%%%%%%%%%
        [Xm, Ym] = meshgrid(linspace(0,14000,100),linspace(0,-14000,100)  );
        fix = scatteredInterpolant(MESH.NODES(1,:)',MESH.NODES(2,:)', U(1:2:end-1));
        fiy = scatteredInterpolant(MESH.NODES(1,:)',MESH.NODES(2,:)', U(2:2:end) );
        
        Xmi = fix(Xm,Ym);
        Ymi = fiy(Xm,Ym);
        
        figure(2),clf
        toplot = sqrt( U(1:2:end-1).^2 + U(2:2:end).^2 ) ;
        toplot = toplot(MESH.ELEMS(1:3,:));
        trisurf(reshape(1:3*nel,3, nel)', MESH.NODES(1,MESH.ELEMS(1:3,:)), MESH.NODES(2,MESH.ELEMS(1:3,:)), zeros(size(MESH.NODES(1,MESH.ELEMS(1:3,:)))), toplot);
        view(2)
        shading interp
        
        hold on
        rad = 2000;
        theta    = linspace(pi/2,-pi/2,40);
        xx_incl = rad * cos(theta);
        yy_incl = rad * sin(theta) - 7000;
        xx_incl([1 end]) = xx_incl([1 end]) + 1;
        hlines = streamline(Xm, Ym, Xmi, Ymi,xx_incl,yy_incl);
        set(hlines,'LineWidth',1,'Color','k')
        
        Xm3 = zeros(size(Xm,1),size(Xm,2),2); Ym3 = Xm3; Zm3 = Xm3;
        Xm3(:,:,1) = Xm; Xm3(:,:,2) = Xm;
        Ym3(:,:,1) = Ym; Ym3(:,:,2) = Ym;
        Zm3(:,:,2) = Zm3(:,:,1)+1;
        
        Xm3i = zeros(size(Xmi,1),size(Xmi,2),2); Ym3i = Xm3i; Zm3i = Xm3i;
        Xm3i(:,:,1) = Xmi; Xm3i(:,:,2) = Xmi;
        Ym3i(:,:,1) = Ymi; Ym3i(:,:,2) = Ymi;
        Zm3i(:,:,2) = Zm3i(:,:,1)+1;
        
        Ampv = sqrt(Xm3i.^2 + Ym3i.^2);
        Xm3i = Xm3i./Ampv;
        Ym3i = Ym3i./Ampv;
        
        xcone = [];
        ycone = [];
        for il = 1:length(hlines)
            
            tmp = hlines(il);
            xy  = [tmp.XData(:) tmp.YData(:)] ;
            xcone = [xcone xy([1 250],1)' xy(end,1)];
            ycone = [ycone xy([1 250],2)' xy(end,2)];
        end
        coneplot(Xm3,Ym3,Zm3,  Xm3i, Ym3i, Zm3i ,xcone, ycone, 0*xcone,2)
        axis equal, box on, grid on
        xlim([-500 10000])
        ylim([-10000 500])
        
        c = colorbar;
        c.Label.String = 'Displacement amplitude [m]';
        c.Label.FontSize = 10;
        c.Label.FontWeight = 'bold';
        caxis([.5 3.5])
        xlabel('Distance [m]','FontSize',10,'FontWeight','bold')
        ylabel('Depth [m]' ,'FontSize',10,'FontWeight','bold')
        title('Numerical displacement field','FontSize',10,'FontWeight','bold')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        figure(3), clf, hold on
        % horizontal surface disp
        surf_nod = unique(MESH.SEGMENTS(:, find(MESH.segment_markers==3)));
        UX = U(1:2:end);
        plot(MESH.NODES(1,surf_nod),UX(surf_nod),'kd')
        % vertical surface disp
        surf_nod = unique(MESH.SEGMENTS(:, find(MESH.segment_markers==3)));
        UY = U(2:2:end);
        plot(MESH.NODES(1,surf_nod),UY(surf_nod),'ks','MarkerFaceColor',[0.5,0.5,0.5])
        grid on, box on
        title('Surface displacements','FontSize',10,'FontWeight','bold')
        xlabel('Distance [m]','FontSize',10,'FontWeight','bold')
        ylabel('Displacement [m]' ,'FontSize',10,'FontWeight','bold')
        xlim([0 10000])
        drawnow
        
        
        %%% TESTING ELASTICITY SURF DEFORMATION FROM EQ.2 GERBAULT 2002 CORRECTED VERSION !!!
        
        H = MESHARG.Mag_depth;
        R = MESHARG.Mag_rad;
        RoverH_val = R/H;
        fun2solv = @(alpha) RoverH_val - ( (2*alpha)/(1+alpha^2) );
        alpha = fzero(fun2solv, 0);
        
        gamma  = (1-alpha^2)/(1+alpha^2);
        mu     = PARAM.ep_G;
        dP     = BC.i_press;
        x      = linspace(0,10000 ,5000);
        
        Ux_ana = (3*dP/2/mu) * (x.*R^2) ./ (x.^2 + gamma^2*H^2);
        Uy_ana = (3*dP/2/mu) * (H*R^2) ./ (x.^2 + gamma^2*H^2);
        
        figure(3), hold on
        plot(x, Ux_ana,'k-','LineWidth',2)
        plot(x, Uy_ana,'k-','LineWidth',2)
        grid on, box on
        legend('Horizontal surf. disp.','Vertical surf. disp.','Analytical solutions')
        
        % Saving errors
        y  = 0*x; 
        Ux_fem = fix(x, y);
        Uy_fem = fiy(x, y);
        erx(ii) = sum(abs(abs(Ux_fem) - abs(Ux_ana) )) / sum(abs(Ux_ana)) *100;
        ery(ii) = sum(abs(abs(Uy_fem) - abs(Uy_ana) )) / sum(abs(Uy_ana)) *100 ;
        
    end
    
end


figure(4), clf
plot(Lres*100, erx,'-kd')
grid on, box on
set(gca,'XScale','log','YScale','log')
hold on
plot(Lres*100, ery,'-ks', 'MarkerFaceColor',[.5 .5 .5])
grid on, box on
set(gca,'XScale','log','YScale','log')
title('Error estimation','FontSize',10,'FontWeight','bold')
ylabel('Numerical error [%]','FontSize',10,'FontWeight','bold')
xlabel('Lateral and vertical limits of the model [km]' ,'FontSize',10,'FontWeight','bold')
legend('Horizontal surf. disp. error','Vertical surf. disp. error')
drawnow





