%% SCRIPT TO EXECUTE 2D ELASTO-PLASTIC BENCHMARK OF
%% HILL's THICK CYLINDER COLLAPSE:

% Part of Open-GeoNabla, copyright GPLv3, 2018
% https://github.com/albansouche/Open-GeoNabla/
% Physics of Geological Processes (PGP) , The NJORD Centre, Dept of Geosciences, University of Oslo
% Author: Alban Souche

% SHORT DETAILS ON SETUP AND SOLVER
% - Cylinder of a given thickness
% - Elasto-platic solver with Von Mises Elasto-plastic constitutive model

% Solution from:
% Hill, Rodney. The mathematical theory of plasticity. Oxford: The Clarendon Press, 1950
% also in (de Souza Neto, E. A., Peric, D., & Owen, D. R., 2011: Computational methods for plasticity: theory and applications. John Wiley & Sons.)


%% CLEARING AND INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(0, 'DefaultFigureRenderer', 'zbuffer');
clear

restoredefaultpath
path_ext_libr = '../../../CODES/ext_libr';
path_2D_libr  = '../../../CODES/2D_libr';
addpath(genpath(path_ext_libr))
addpath(genpath(path_2D_libr))

%% MESH GENERATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dimension
ndim                 = 2;
% input arguments
mesharg.radius_in    = 0.1;  % int radius
mesharg.radius_ext   = 0.2;  % ext radius
mesharg.angle_sym    = 90*pi/180;
mesharg.max_tri_area = 1e-6; % maximum triangle area
mesharg.type_el      = 'tri7'; % 'tri3' or 'tri7'
% generate mesh
MESH = generate_mesh_cylinder(mesharg);
% reordering
[MESH,~,~] = mesh_reorder_amd(MESH);
% model parameter
nnod         = size(MESH.NODES,2);
[nnodel,nel] = size(MESH.ELEMS);
fprintf(1, ['\n Number of nodes:   ', num2str(nnod)]);
fprintf(1, ['\n Number of elems:   ', num2str(nel),'\n']);

% plotting meshes
figure(1); clf;
toplot = repmat(MESH.elem_markers,3,1);
trisurf(reshape(1:3*nel,3, nel)', MESH.NODES(1,MESH.ELEMS(1:3,:)), MESH.NODES(2,MESH.ELEMS(1:3,:)), zeros(size(MESH.NODES(1,MESH.ELEMS(1:3,:)))), toplot);
title('mesh'); view(2); axis image; axis off; axis equal

%% Solver
SOLVER.nelblo      = 10000; % number of ele per block for matrix assemblage (MILAMIN approach)
SOLVER.nb_cpu      = 2;     % number of cpu to use for parallel calculations
SOLVER.solver_type = 'direct_sym_chol'; % 'direct_sym_chol' / 'direct_sym' / 'direct_nonsym' / (develpment: 'pcg'/'AGMG')

%% MESH manipulation: segments and normals
MESH  = seg_manip(MESH, mesharg.type_el, SOLVER);
MESH  = NormalEdges(MESH);

%% Integration and evaluation points
if nnodel==3
    SOLVER.nip  = 1;
    SOLVER.neva = 1; % Constant discontinuous strain field
elseif nnodel==7
    SOLVER.nip  = 6;
    SOLVER.neva = 6; % Linear discontinuous strain field
end
[SOLVER.ieuv, ~] = ip_triangle(SOLVER.neva);


%% PHYSICS Elasto Plastic domain %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PARAM.ep_S_y    = 0.24 * 1e9; % Yield stress [Pa]
PARAM.ep_E      = 210  * 1e9; % Young's modulus [Pa]
PARAM.ep_nu     = 0.499; % Poisson's ratio
PARAM.ep_lambda = PARAM.ep_E.*PARAM.ep_nu./((1+PARAM.ep_nu).*(1-2*PARAM.ep_nu));
PARAM.ep_G      = PARAM.ep_E./(2*(1+PARAM.ep_nu)); % shear modulus
PARAM.ep_K      = PARAM.ep_E./(3*(1-2*PARAM.ep_nu)); % bulk modulus
PARAM.ep_H      = 0; % isotropic hardening
PARAM.ep_rho    = 0;      % density [kg/m3] % set to 0 for VM material
PARAM.g         = [0; 0]; % gravity [m/s2]  % set to 0 for VM material
PARAM.ep_phi    = 0;  % set to 0 for VM material
PARAM.ep_xsi    = 0;  % set to 0 for VM material


%% Boundary conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bc_ind_ux = unique(MESH.SEGMENTS(:,MESH.segment_markers==2));
bc_val_ux = zeros(size(bc_ind_ux));
bc_ind_uy = unique(MESH.SEGMENTS(:,MESH.segment_markers==4));
bc_val_uy = zeros(size(bc_ind_uy));
BC.ind    = [ndim*(bc_ind_ux(:)-1)+1; ndim*(bc_ind_uy(:)-1)+2];
BC.val    = [bc_val_ux(:); bc_val_uy(:)];
BC2       = BC;        % for Newton-Raphson correction calculation
BC2.val   = 0*BC2.val; % for Newton-Raphson correction calculation
hold on
plot(MESH.NODES(1,bc_ind_ux), MESH.NODES(2,bc_ind_ux), 'yo')
plot(MESH.NODES(1,bc_ind_uy), MESH.NODES(2,bc_ind_uy), 'bo')


%% Analytical solution (Hill's) FROM De Suza book (Tresca with modified Yield stress to approxiamte von Mises) %%%%%%%%
a = mesharg.radius_in;
b = mesharg.radius_ext;
Y = 2 * PARAM.ep_S_y / (sqrt(3));
Po    = Y * (1-(a^2/b^2)) / 2;
Pinf  = Y * log(b/a);
Pload_e = linspace(0,Pinf,100);
Pload_e(end) = [];
ub_e  = 2*Pload_e*b*(1-PARAM.ep_nu^2)/PARAM.ep_E/((b^2/a^2)-1);
c = linspace(a,b,1000);
Pload_p = Y * (log(c/a) + 0.5*(1-(c.^2/b^2)) );
ub_p = Y*(c.^2)*(1-PARAM.ep_nu^2)/PARAM.ep_E/b;
figure(2); clf
hold on, plot(ub_e,Pload_e,'k'), plot(ub_p,Pload_p,'r')
% Node for Numeric vs Analytical comparison
indx_ext_rad = bc_ind_uy(MESH.NODES(1,bc_ind_uy)==max(MESH.NODES(1,bc_ind_uy)));

%% Internal pressure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_max    = 0.99 * Pinf; %[Pa]
P_ini    = 0;
nb_incre = 20;
P_incre  = (P_max-P_ini)/nb_incre;

%% dU elastic increment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Elastic tangent matix
D_el_ip = Tangent_Elastoplastic_VonMises(MESH, [], [], [], PARAM, SOLVER.nip);
% Traction BC (Fext)
SIGMA = zeros(3,nnod);
seg_id   = 1;
tmp_nods = unique(MESH.SEGMENTS(:,MESH.segment_markers==seg_id));
SIGMA(1,tmp_nods) = - P_incre;
SIGMA(2,tmp_nods) = - P_incre;
fext  = ContourTraction(MESH, SIGMA, seg_id, []); % Rhs Traction vector (external forces)
clear SIGMA
figure(1),hold on
quiver(MESH.NODES(1,:),MESH.NODES(2,:),fext(1:2:end-1)',fext(2:2:end)')
plot(MESH.NODES(1,tmp_nods),MESH.NODES(2,tmp_nods),'dr' )
axis equal
drawnow
% Assemble stiffness matrix and rhs
elast = 1; % '1'== elastic only ; '0'==elastoplastic
[A_all, ~, InvJ_all, DetJ_all] = matrix_entries_EP(MESH, D_el_ip, PARAM, SOLVER, elast);
%Solve for elastic predicator increment
dU_e    = mechanical_solver(MESH, A_all, BC, fext, SOLVER);
dE_e    = strain_tot_calc(MESH, InvJ_all, dU_e, SOLVER.ieuv, SOLVER.nelblo);
dE_e.xy = 2*dE_e.xy; % gamma convention

%% Initialisation
U         = zeros(ndim*nnod,1); % total displacement
Etot.xx   = zeros(SOLVER.neva, nel);   % total strain
Etot.yy   = zeros(SOLVER.neva, nel);
Etot.zz   = zeros(SOLVER.neva, nel);
Etot.xy   = zeros(SOLVER.neva, nel);
Epl.bar_n = zeros(SOLVER.neva,nel);    % plastic strain
Epl.xx_n  = zeros(SOLVER.neva,nel);
Epl.yy_n  = zeros(SOLVER.neva,nel);
Epl.zz_n  = zeros(SOLVER.neva,nel);
Epl.xy_n  = zeros(SOLVER.neva,nel);
figure(4); clf
figure(5); clf


%% INCREMENTAL LOAD PROCEDURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_it_nr   = 100;
rel_res_tol = 1e-9;

for iincre = 1:nb_incre

    %%% External load increment  (used for NR residual)%%%%%%%%%%%%%%%%%%%%
    Fext = iincre * fext;

    %%% Initialization
    rel_res   	= NaN(1, max_it_nr+1);
    rel_res(1) 	= 1;
    d_gamma = zeros(size(Etot.xx));
    Etrial  = Etot;
    dU      = dU_e;
    Utrial  = U + dU_e;

    %%% Iterations
    for i = 1:max_it_nr

        if i==1
            dEtrial = dE_e;
        else
            dEtrial      = strain_tot_calc(MESH, InvJ_all, dU, SOLVER.ieuv, SOLVER.nelblo);
            dEtrial.xy   = 2*dEtrial.xy; % gamma convention
        end

        % Add strain increment
        Etrial.xx = Etot.xx + dEtrial.xx;
        Etrial.yy = Etot.yy + dEtrial.yy;
        Etrial.zz = Etot.zz + dEtrial.zz;
        Etrial.xy = Etot.xy + dEtrial.xy;
        Etrial.ii = Etrial.xx+Etrial.yy+Etrial.zz;
        % Calculate Pressure
        P       = PARAM.ep_K*Etrial.ii;
        % Calculate deviatoric stress
        Sdev.xx = 2*PARAM.ep_G * ( Etrial.xx - Etrial.ii/3 );
        Sdev.yy = 2*PARAM.ep_G * ( Etrial.yy - Etrial.ii/3 );
        Sdev.zz = 2*PARAM.ep_G * ( Etrial.zz - Etrial.ii/3 );
        Sdev.xy = 2*PARAM.ep_G * ( Etrial.xy/2 );

        % Calculate J2
        J2 = (Sdev.xx.^2 + Sdev.yy.^2 + Sdev.zz.^2)/2 + Sdev.xy.^2;
        % Calculate Von Mises effective stress
        Q  = sqrt(3*J2);

        % Check for plasticity consistency %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        F      = Q - PARAM.ep_S_y;
        indx_p = F>0;

        d_gamma = zeros(SOLVER.nip,nel);

        if find(indx_p)

            % Return mapping
            d_g = Return_mapping_VonMises(MESH, F, PARAM, indx_p);
            d_gamma(indx_p) = d_g;

            % Consistent tangent modulus (using non-corrected stress tensor!)
            D_el_ip = Tangent_Elastoplastic_VonMises(MESH, indx_p, d_gamma, Sdev, PARAM, SOLVER.nip);

            % Update stiffness matrix
            elast = 0;
            [A_all, ~, ~, ~] = matrix_entries_EP(MESH, D_el_ip, PARAM, SOLVER, elast);

            % Apply correction to the stress tensor
            fact = 1 - 3*PARAM.ep_G*d_gamma./Q;
            Sdev.xx = fact.*Sdev.xx;
            Sdev.yy = fact.*Sdev.yy;
            Sdev.zz = fact.*Sdev.zz;
            Sdev.xy = fact.*Sdev.xy;

            % Compute residual
            Fint = internal_force(MESH, InvJ_all, DetJ_all, Sdev, P, SOLVER);
            res          = Fext - Fint;
            res(BC.ind)  = 0;
            rel_res(i+1) = norm(res)/norm(Fext);
            fprintf(1, ['Newton-Raphson iteration: ', num2str(i,'%3.0d'), '    relative residual: ', num2str(rel_res(i+1),'%1.2e'),'\n']);

            if rel_res(i+1)<rel_res_tol
                fprintf(1, '\n');
                break
            else
                ddU = mechanical_solver(MESH, A_all, BC2, res, SOLVER);
                dU  = dU + ddU; % Newton-Raphson correction
            end

        else
            break
        end

        Utrial = U + dU; % Update trial solution

        Ux = Utrial(1:2:end-1);
        Uy = Utrial(2:2:end);
        Uamp = sqrt(Ux.^2 + Uy.^2);
        figure(2); hold on;
        plot(Uamp(indx_ext_rad),iincre*P_incre,'xk')
        drawnow

    end

    % Update converged solution
    Etot = Etrial;
    U = U + dU;

    % Plotting
    Ux = U(1:2:end-1);
    Uy = U(2:2:end);
    Uamp = sqrt(Ux.^2 + Uy.^2);
    figure(2); hold on;
    plot(Uamp(indx_ext_rad),iincre*P_incre,'sk')
    xlabel('External radius displacement [m]','FontSize',10,'FontWeight','bold')
    ylabel('Inner pressure load [Pa]','FontSize',10,'FontWeight','bold')
    grid on, box on
    drawnow

    if ismember(iincre, [10 15 20])

        figure(4), hold on
        F = scatteredInterpolant(MESH.NODES(1,:)',MESH.NODES(2,:)', Ux);
        a = 0.1;
        b = 0.2;
        press     = iincre*P_incre;
        rad       = linspace(a, b, 20);
        u_num_rad = F(rad, 0*rad);
        plot(rad, u_num_rad,'kd','MarkerFaceColor',[0.5,0.5,0.5])
        ez        = 0*((1-(2*PARAM.ep_nu))*press) / (PARAM.ep_E * ((b*b)/(a*a) - 1) );
        ur_e_Hill = - PARAM.ep_nu.*rad*ez ...
            + ((1+PARAM.ep_nu)*press) / (PARAM.ep_E * ( (b*b)/(a*a)-1) ) * ( (1-2*PARAM.ep_nu).*rad + ((b*b)./rad) );
        Y   = 2*PARAM.ep_S_y / (sqrt(3));  % vonMises approx. to Tresca solution
        fun = @(c) ( - press + (Y*log(c/a)+Y/2*(1-(c*c)/(b*b))) ) ;
        c   = fzero(fun, a);
        if c>a
            ur_ep_Hill = - PARAM.ep_nu.*rad*ez ...
                + ((1+PARAM.ep_nu)*Y*c*c) / (2*PARAM.ep_E*b*b) * ( (1-2*PARAM.ep_nu).*rad + (b*b./rad) );
            plot(rad, ur_ep_Hill,'-k')
        else
            plot(rad, ur_e_Hill,'-k')
        end
        xlabel('Radius [m]','FontSize',10,'FontWeight','bold')
        ylabel('Radial displacement [m]','FontSize',10,'FontWeight','bold')
        grid on, box on
        drawnow
    end

    figure(3); clf;
    if nnodel==3
        toplot = repmat(P,3,1);
    else
        toplot = repmat(mean(P),3,1);
    end
    trisurf(reshape(1:3*nel,3, nel)', MESH.NODES(1,MESH.ELEMS(1:3,:)), MESH.NODES(2,MESH.ELEMS(1:3,:)), zeros(size(MESH.NODES(1,MESH.ELEMS(1:3,:)))), toplot);
    title('P'); view(2); colorbar; axis image; axis off
    caxis([min(toplot(:)) max(toplot(:))]);
    shading interp

    if find(indx_p)
        figure(5); hold on
        loglog(1:max_it_nr,rel_res(2:end),'-dk','MarkerFaceColor',[0.5,0.5,0.5])
        xlabel('Number of iteration','FontSize',10,'FontWeight','bold')
        ylabel('Relative residual','FontSize',10,'FontWeight','bold')
        set(gca,'Xscale','linear','Yscale','log')
        grid on, box on
        drawnow
    end

end
