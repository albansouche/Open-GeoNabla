function [U, dP_c_i, dP_c, dP_cm] = EP_DruckerPrager_Critic_dP(MESH, PARAM, BC, P_litho, Fext, SOLVER, POST)

%   Implementation of Compressible Elasto Plasticity
%   Code build upon MILAMIN Version 1.0.1, Mutils 0.4.2, HYPLAS_v2.0, FOLDER 1.0
%   By ALBAN SOUCHE, Physics of Geological Processes, Njord Center, Departement of Geosciences, University of Oslo, January 2017

%% Mesh parameters
[ndim, nnod] = size(MESH.NODES);
[~   ,  nel] = size(MESH.ELEMS);

%% Integration and evaluation points
nip       = SOLVER.nip;
neva      = SOLVER.neva;
[ieuv, ~] = ip_tetra(neva);

%% Boundary conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BC_NR.ind = BC.ep_ind;   % for Newton-Raphson correction calculation
BC_NR.val = 0*BC.ep_val;
BC_e.ind  = BC.ep_ind;   % for elastic predicator increment
BC_e.val  = BC.ep_val;

%% dU elastic increment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Elastic tangent matix
indx_p  = false(nip,nel);
D_el_ip = Tangent_Elastoplastic_DruckerPrager(MESH, indx_p, [], [], [], PARAM, SOLVER.nip);

% Traction BC (Fext)
nb_incre = SOLVER.ep_nb_incre;
fext     = Fext/nb_incre; % incremental external load vectorSIGMA = zeros(3,nnod);
% Assemble stiffness matrix and rhs
elast = 1; % '1'== elastic only ; '0'==elastoplastic
[A_all, InvJ_all, DetJ_all] = matrix_entries_EP(MESH, D_el_ip, PARAM, SOLVER, elast);

%Solve for elastic predicator increment
SOLVER_e= SOLVER;
%SOLVER_e.solver_type = 'direct_sym_chol';
SOLVER_e.solver_type = SOLVER.solver_type;
dU_e    = mechanical_solver(MESH, A_all, BC_e, fext, SOLVER_e); % fext increment as rhs
dE_e    = strain_tot_calc(MESH, InvJ_all, dU_e, ieuv, SOLVER.nelblo);
dE_e.xy = 2*dE_e.xy; % gamma convention
dE_e.yz = 2*dE_e.yz; % gamma convention
dE_e.zx = 2*dE_e.zx; % gamma convention

% %% Initialisation
U       = zeros(ndim*nnod,1); % total displacement
Etot.xx = zeros(neva, nel);   % total strain
Etot.yy = zeros(neva, nel);
Etot.zz = zeros(neva, nel);
Etot.xy = zeros(neva, nel);
Etot.yz = zeros(neva, nel);
Etot.zx = zeros(neva, nel);
Epl.bar = zeros(neva,nel); % accumulated plastic strain (scalar)
Epl.xx  = zeros(neva,nel); % plastic strain tensor
Epl.yy  = zeros(neva,nel);
Epl.zz  = zeros(neva,nel);
Epl.xy  = zeros(neva,nel);
Epl.yz  = zeros(neva,nel);
Epl.zx  = zeros(neva,nel);
nelblo  = SOLVER.nelblo;
rel_res_tol = SOLVER.ep_rel_res_tol;


%% INCREMENTAL LOAD PROCEDURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_it_nr = SOLVER.ep_max_NR_it;
iincre = 1;
out_loop = 0;

% Checking if failure is at interface or not
indx7 = unique(find(MESH.node_markers==7));
elem7 = find(sum(ismember(MESH.ELEMS,indx7)));
dP_c_i=0;

while and(iincre <= nb_incre, out_loop==0)
    
    %%% External load increment  (used for NR residual)%%%%%%%%%%%%%%%%%%%%
    Fext = iincre*fext;
    fprintf(1, ['Load increment ',num2str(iincre),' MPa / ',num2str(nb_incre),' MPa \n']);
    
    %%% Initialization
    rel_res   	= NaN(1, max_it_nr+1);
    rel_res(1) 	= 1;
    dU = dU_e;
    
    %%% Iterations
    for i = 1:max_it_nr
        
        TimeIncr = tic;
        
        if i==1
            dEtrial = dE_e;
        else
            dEtrial      = strain_tot_calc(MESH, InvJ_all, dU, ieuv, nelblo);
            dEtrial.xy   = 2*dEtrial.xy; % gamma convention
            dEtrial.yz   = 2*dEtrial.yz; % gamma convention
            dEtrial.zx   = 2*dEtrial.zx; % gamma convention
        end
        
        % Trial strain
        Etrial.xx = Etot.xx + dEtrial.xx;
        Etrial.yy = Etot.yy + dEtrial.yy;
        Etrial.zz = Etot.zz + dEtrial.zz;
        Etrial.xy = Etot.xy + dEtrial.xy;
        Etrial.yz = Etot.yz + dEtrial.yz;
        Etrial.zx = Etot.zx + dEtrial.zx;
        Etrial.ii = Etrial.xx+Etrial.yy+Etrial.zz;
        % Trail accumulated plastic strain
        Epl.bar_trial = Epl.bar;
        % Trial presure
        P = bsxfun(@times,PARAM.ep_K(MESH.elem_markers),Etrial.ii) - P_litho;  %% add (-P_litho) negative because fo convention
        
        % Calculate deviatoric stress
        Sdev.xx = bsxfun(@times, 2*PARAM.ep_G(MESH.elem_markers) , ( Etrial.xx - Etrial.ii/3 ) );
        Sdev.yy = bsxfun(@times, 2*PARAM.ep_G(MESH.elem_markers) , ( Etrial.yy - Etrial.ii/3 ) );
        Sdev.zz = bsxfun(@times, 2*PARAM.ep_G(MESH.elem_markers) , ( Etrial.zz - Etrial.ii/3 ) );
        Sdev.xy = bsxfun(@times, 2*PARAM.ep_G(MESH.elem_markers) , ( Etrial.xy/2 ) );
        Sdev.yz = bsxfun(@times, 2*PARAM.ep_G(MESH.elem_markers) , ( Etrial.yz/2 ) );
        Sdev.zx = bsxfun(@times, 2*PARAM.ep_G(MESH.elem_markers) , ( Etrial.zx/2 ) );
        
        % Trial J2
        J2 = (Sdev.xx.^2 + Sdev.yy.^2 + Sdev.zz.^2)/2 + Sdev.xy.^2 + Sdev.yz.^2 + Sdev.zx.^2;
        
        % Check for plasticity consistency %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cohe   = bsxfun(@plus, PARAM.ep_Co(MESH.elem_markers), PARAM.ep_H*Epl.bar);
        F      = sqrt(J2) + PARAM.ep_eta*P - PARAM.ep_xi*cohe;
        indx_p = F>0;
        d_gamma = zeros(nip,nel);
            
        if find(indx_p)
            
            % Checking if failure is at interface or not
            if find(ismember(find(sum(indx_p,1)),elem7))
                figure(1),hold on
                tetramesh(MESH.ELEMS(1:4, find(sum(indx_p,1)) )' , MESH.NODES'  )
                fprintf(1,'failure at interface \n \n \n \n');
                out_loop = 1;
                dP_c = iincre;
                tmp = -(P + P_litho);
                dP_cm = max(tmp(:)) ;
                break
            else 
                figure(1),hold on
                tetramesh(MESH.ELEMS(1:4, find(sum(indx_p,1)) )' , MESH.NODES'  )
                fprintf(1,'failure not at interface \n \n \n \n');
                if dP_c_i==0
                    dP_c_i = iincre;
                end    
            end
            
            % Return mapping with update of deviatoric stress tensor and P
            [d_g, Sdev, P, Epl, indx_apex] = Return_mapping_DruckerPrager(MESH, F, PARAM, indx_p, J2, Sdev, P, Epl);
            d_gamma(indx_p) = d_g;
            
            % Consistent tangent modulus (using non-corrected strain tensor!)
            D_el_ip = Tangent_Elastoplastic_DruckerPrager(MESH, indx_p, indx_apex, d_gamma, Etrial, PARAM, nip);
            
            % Update stiffness matrix
            elasto = 0;
            [A_all, ~, ~] = matrix_entries_EP(MESH, D_el_ip, PARAM, SOLVER, elasto);
            
            % Compute residual
            Fint         = internal_force(MESH, InvJ_all, DetJ_all, Sdev, P + P_litho, SOLVER);
            res          = Fext - Fint;
            res(BC_NR.ind)  = 0;
            rel_res(i+1) = norm(res)/norm(Fext);
            
            if rel_res(i+1)<rel_res_tol
                % print info to screen
                fprintf(1, ['Load increment (',num2str(iincre),'/',num2str(nb_incre),')    NR iterations: ', num2str(i,'%3.0d'), '    relative residual: ', num2str(rel_res(i+1),'%1.2e'),'    time(s): ',num2str(toc(TimeIncr)) ,'\n\n']);
                break
            else
                % Calculate NR increment
                ddU = mechanical_solver(MESH, A_all, BC_NR, res, SOLVER); % res as rhs
                alpha = 1; % Full Newton-Raphson update
                dU  = dU + alpha*ddU; % Updated solution
            end
            
            % print info to screen
            fprintf(1, ['Load increment (',num2str(iincre),'/',num2str(nb_incre),')    NR iterations: ', num2str(i,'%3.0d'), '    relative residual: ', num2str(rel_res(i+1),'%1.2e'),'    time(s): ',num2str(toc(TimeIncr)) ,'\n']);

        else
            break
        end
        
    end
    
    iincre = iincre+1;
    
    % Update converged solution
    Etot    = Etrial;
    Epl.bar = Epl.bar_trial;
    U       = U + dU;
    
    %%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if or(mod(iincre, POST.Hz_action)==0, iincre==1)
        
        type_el = 'tet4';
        elems_tmp = MESH.ELEMS(1:4,:);
        
        if POST.action(2)==1
            % Trial J2
            J2 = (Sdev.xx.^2 + Sdev.yy.^2 + Sdev.zz.^2)/2 + Sdev.xy.^2 + Sdev.yz.^2 + Sdev.zx.^2;
            cohe   = bsxfun(@plus, PARAM.ep_Co(MESH.elem_markers), PARAM.ep_H*Epl.bar);
            F      = sqrt(J2) + PARAM.ep_eta*P - PARAM.ep_xi*cohe;
            if SOLVER.neva==1
                tmp = repmat(F,4,1);
                tmp_n = accumarray(elems_tmp(:),tmp(:));
                tmp_n = tmp_n./MESH.EperN;
            elseif SOLVER.neva==4
                tmp_n = accumarray(elems_tmp(:),F(:));
                tmp_n = tmp_n./MESH.EperN;
            end
            paraview_write(type_el ,MESH.NODES, elems_tmp, [], tmp_n', [POST.folder_name,'/Yield_fun'], iincre)
        end
        
        if POST.action(3)==1
            if SOLVER.neva==1
                tmp = repmat(P+P_litho,4,1);
                tmp_n = accumarray(elems_tmp(:),tmp(:));
                tmp_n = tmp_n./MESH.EperN;
            elseif SOLVER.neva==4
                tmp_n = accumarray(elems_tmp(:),P(:)+P_litho(:));
                tmp_n = tmp_n./MESH.EperN;
            end
            paraview_write(type_el ,MESH.NODES, elems_tmp, [], tmp_n', [POST.folder_name,'/Pressure'], iincre)
        end
        
        if POST.action(4)==1
            if SOLVER.neva==1
                tmp = repmat(J2,4,1);
                tmp_n = accumarray(elems_tmp(:),tmp(:));
                tmp_n = tmp_n./MESH.EperN;
            elseif SOLVER.neva==4
                tmp_n = accumarray(elems_tmp(:),J2(:));
                tmp_n = tmp_n./MESH.EperN;
            end
            paraview_write(type_el ,MESH.NODES, elems_tmp,  [], tmp_n', [POST.folder_name,'/J2'], iincre)
        end
        
        if POST.action(5)==1
            E2 = ( (Etot.xx- Etrial.ii/3).^2 + (Etot.yy- Etrial.ii/3).^2 + (Etot.zz- Etrial.ii/3).^2 )/2 + (Etrial.xy/2).^2  + (Etrial.yz/2).^2  + (Etrial.zx/2).^2 ;
            if SOLVER.neva==1
                tmp = repmat(E2,4,1);
                tmp_n = accumarray(elems_tmp(:),tmp(:));
                tmp_n = tmp_n./MESH.EperN;
            elseif SOLVER.neva==4
                tmp_n = accumarray(elems_tmp(:),E2(:));
                tmp_n = tmp_n./MESH.EperN;
            end
            paraview_write(type_el ,MESH.NODES, elems_tmp,  [], tmp_n', [POST.folder_name,'/Strain_E2'], iincre)
        end
        
        if POST.action(6,1)==1
            if SOLVER.neva==1
                tmp = repmat(Etrial.ii,4,1);
                tmp_n = accumarray(elems_tmp(:),tmp(:));
                tmp_n = tmp_n./MESH.EperN;
            elseif SOLVER.neva==4
                tmp_n = accumarray(elems_tmp(:),Etrial.ii(:));
                tmp_n = tmp_n./MESH.EperN;
            end
            paraview_write(type_el ,MESH.NODES, elems_tmp,  [], tmp_n', [POST.folder_name,'/Strain_vol'], iincre)
        end
        
        if POST.action(7)==1
            paraview_write(type_el ,MESH.NODES, elems_tmp, MESH.elem_markers, [U(1:3:end-2)'; U(2:3:end-1)'; U(3:3:end)'] , [POST.folder_name,'/U_vec'], iincre)
        end

        
    end
    
end






