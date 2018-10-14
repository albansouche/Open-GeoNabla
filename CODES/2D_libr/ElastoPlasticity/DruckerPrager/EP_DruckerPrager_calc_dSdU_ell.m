function U = EP_DruckerPrager_calc_dSdU_ell(MESH, PARAM, BC, P_litho, Fext, SOLVER, POST)

%   Implementation of Compressible Elasto Plasticity
%   Code build upon MILAMIN Version 1.0.1, Mutils 0.4.2, FOLDER 1.0
%   By ALBAN SOUCHE, PGP UiO, August 2016

%% Mesh parameters
[ndim, nnod] = size(MESH.NODES);
[nnodel,nel] = size(MESH.ELEMS);

%% Integration and evaluation points
nip       = SOLVER.nip;
neva      = SOLVER.neva;
[ieuv, ~] = ip_triangle(neva);

%% Boundary conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BC_NR.ind = BC.ep_ind;   % for Newton-Raphson correction calculation
BC_NR.val = 0*BC.ep_val;
BC_e.ind = BC.ep_ind;   % for elastic predicator increment
BC_e.val = BC.ep_val;

%% dU elastic increment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Elastic tangent matix
indx_p  = false(nip,nel);
D_el_ip = Tangent_Elastoplastic_DruckerPrager(MESH, indx_p, [], [], [], PARAM, SOLVER.nip);

% Traction BC (Fext)
nb_incre = SOLVER.ep_nb_incre;
fext     = Fext/nb_incre; % incremental external load vector SIGMA = zeros(3,nnod);
% Assemble stiffness matrix and rhs
elast = 1; % '1'== elastic only ; '0'==elastoplastic
[A_all, ~, InvJ_all, DetJ_all] = matrix_entries_EP(MESH, D_el_ip, PARAM, SOLVER, elast);

%Solve for elastic predicator increment
SOLVER_e= SOLVER;
SOLVER_e.solver_type = 'direct_sym_chol';
dU_e    = mechanical_solver(MESH, A_all, BC_e, fext, SOLVER_e); % fext increment as rhs
dE_e    = strain_tot_calc(MESH, InvJ_all, dU_e, ieuv, SOLVER.nelblo);
dE_e.xy = 2*dE_e.xy; % gamma convention

% %% Initialisation
U       = zeros(ndim*nnod,1); % total displacement
Etot.xx = zeros(neva, nel);   % total strain
Etot.yy = zeros(neva, nel);
Etot.zz = zeros(neva, nel);
Etot.xy = zeros(neva, nel);
Epl.bar = zeros(neva,nel); % accumulated plastic strain (scalar)
Epl.xx  = zeros(neva,nel); % plastic strain tensor
Epl.yy  = zeros(neva,nel);
Epl.zz  = zeros(neva,nel);
Epl.xy  = zeros(neva,nel);
nelblo  = SOLVER.nelblo;
rel_res_tol = SOLVER.ep_rel_res_tol;

%% INCREMENTAL LOAD PROCEDURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_it_nr = SOLVER.ep_max_NR_it;
dUmod     = 0;

for iincre = 1:nb_incre
    
    %%% External load increment  (used for NR residual)%%%%%%%%%%%%%%%%%%%%
    Fext = iincre*fext;
    fprintf(1, ['Increment (',num2str(iincre),'/',num2str(nb_incre),')\n']);
    
    %%% Initialization
    rel_res   	= NaN(1, max_it_nr+1);
    rel_res(1) 	= 1;
    dU = dU_e;
        
    %%% Iterations
    for i = 1:max_it_nr
        
        
       if dUmod==0
           
            TimeIncr = tic;
            
            if i==1
                dEtrial = dE_e;
            else
                dEtrial      = strain_tot_calc(MESH, InvJ_all, dU, ieuv, nelblo);
                dEtrial.xy   = 2*dEtrial.xy; % gamma convention
            end
            
            % Trial strain
            Etrial.xx = Etot.xx + dEtrial.xx;
            Etrial.yy = Etot.yy + dEtrial.yy;
            Etrial.zz = Etot.zz + dEtrial.zz;
            Etrial.xy = Etot.xy + dEtrial.xy;
            Etrial.ii = Etrial.xx+Etrial.yy+Etrial.zz;
            % Trail accumulated plastic strain
            Epl.bar_trial = Epl.bar;
            % Trial presure
            P = bsxfun(@times,PARAM.ep_K(MESH.elem_markers),Etrial.ii) - P_litho;  %% add (-P_litho) negative because of convention
            
            % Calculate deviatoric stress
            Sdev.xx = bsxfun(@times, 2*PARAM.ep_G(MESH.elem_markers) , ( Etrial.xx - Etrial.ii/3 ) );
            Sdev.yy = bsxfun(@times, 2*PARAM.ep_G(MESH.elem_markers) , ( Etrial.yy - Etrial.ii/3 ) );
            Sdev.zz = bsxfun(@times, 2*PARAM.ep_G(MESH.elem_markers) , ( Etrial.zz - Etrial.ii/3 ) );
            Sdev.xy = bsxfun(@times, 2*PARAM.ep_G(MESH.elem_markers) , ( Etrial.xy/2 ) );
            
            % Trial J2
            J2 = (Sdev.xx.^2 + Sdev.yy.^2 + Sdev.zz.^2)/2 + Sdev.xy.^2;
            
            % Check for plasticity consistency %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % cohe   = bsxfun(@plus, PARAM.ep_Co(MESH.elem_markers), PARAM.ep_H*Epl.bar);
            cohe   = bsxfun(@plus, PARAM.ep_Co(MESH.elem_markers) + PARAM.noise, PARAM.ep_H*Epl.bar);
            F      = sqrt(J2) + PARAM.ep_eta*P - PARAM.ep_xi*cohe;
            indx_p = F>0; 
            
            d_gamma = zeros(nip,nel);
            
            if find(indx_p)
                     
                
                
                % Return mapping with update of deviatoric stress tensor and P
                [d_g, Sdev, P, Epl, indx_apex] = Return_mapping_DruckerPrager_cohesion_noise(MESH, F, PARAM, indx_p, J2, Sdev, P, Epl);
                d_gamma(indx_p) = d_g;
                
                %%%%%%%%
%                 figure(1),clf, hold on
%                 for iii = 1:length(find(indx_p))
%                     pts_indx = find(indx_p);
%                     pts_indx = pts_indx(iii);
%                     MohrDiag_point(pts_indx,PARAM,Sdev,P,cohe(pts_indx))
%                     drawnow
%                 end
%                 
%                 figure(1),clf, hold on, box on, grid on
%                 J2 = (Sdev.xx.^2 + Sdev.yy.^2 + Sdev.zz.^2)/2 + Sdev.xy.^2;
%                 pts_indx = find(indx_p);
%                 y_line = linspace(0,2*max(sqrt(J2(pts_indx))), 100);
%                 co = PARAM.ep_xi*cohe(pts_indx);
%                 x_line = -(co - y_line)./PARAM.ep_eta;
%                 fact_co = cohe(pts_indx)./PARAM.ep_Co;
%                 plot(x_line, sqrt(3)*y_line .* fact_co ,'-r')
%                 hold on
%                 plot(-P(pts_indx), sqrt(3*J2(pts_indx)) .* fact_co, 'gd')
%                 drawnow
                
%             J2 = (Sdev.xx.^2 + Sdev.yy.^2 + Sdev.zz.^2)/2 + Sdev.xy.^2;
%             F      = sqrt(J2) + PARAM.ep_eta*P - PARAM.ep_xi*cohe;
%             figure(1); clf;
%             if nnodel==3
%                 toplot = repmat(F,3,1);
%             else
%                 toplot = repmat(mean(F),3,1);
%             end
%             trisurf(reshape(1:3*nel,3, nel)', MESH.NODES(1,MESH.ELEMS(1:3,:)), MESH.NODES(2,MESH.ELEMS(1:3,:)), zeros(size(MESH.NODES(1,MESH.ELEMS(1:3,:)))), toplot);
%             title('Yield function F'); view(2); colorbar; axis image; axis off
%             xlim([xminp xmaxp])
%             ylim([yminp ymaxp])
%             shading interp
%             drawnow
                
                %%%%%%%%
                
                % Consistent tangent modulus (using non-corrected strain tensor!)
                D_el_ip = Tangent_Elastoplastic_DruckerPrager(MESH, indx_p, indx_apex, d_gamma, Etrial, PARAM, nip);
                
%                 % Update stiffness matrix
%                 if isequal(PARAM.ep_phi,PARAM.ep_xsi) % associative plasticity rule
                    elast = 0;
                    [A_all, ~, ~, ~] = matrix_entries_EP(MESH, D_el_ip, PARAM, SOLVER, elast);
%                 else
%                     Ael=elasticity_matrix(MESH.NODES,MESH.ELEMS,D_el_ip);
%                     A_all = reshape(Ael,size(Ael,1)*size(Ael,2),size(Ael,3));
%                 end
                
                % Compute residual
                Fint           = internal_force(MESH, InvJ_all, DetJ_all, Sdev, P + P_litho, SOLVER);
                res            = Fext - Fint;
                res(BC_NR.ind) = 0;
                rel_res(i+1)   = norm(res)/norm(Fext);
                
                if rel_res(i+1)<rel_res_tol
                    % print info to screen
                    fprintf(1, ['Load increment (',num2str(iincre),'/',num2str(nb_incre),')    NR iterations: ', num2str(i,'%3.0d'), '    relative residual: ', num2str(rel_res(i+1),'%1.2e'),'    time(s): ',num2str(toc(TimeIncr)) ,'\n\n']);
                    break
                    
                elseif rel_res(i+1)>5e1
                    dUmod=1;    % pass to dU controled solver
                    
                elseif rel_res(i+1)<5e1
                    % Calculate NR increment
                    ddU = mechanical_solver(MESH, A_all, BC_NR, res, SOLVER); % res as rhs
                    if isequal(SOLVER.ep_NR,'NR_line_search') % Line search
                        f = @(s) norm( res_calc_DP(MESH, InvJ_all, DetJ_all, Etot, Epl, PARAM, ieuv, BC_NR, indx_p, dU+s*ddU, Fext, P_litho, SOLVER) );
                        alpha = poly_search(f,norm(res));
                    elseif isequal(SOLVER.ep_NR,'NR_full_update') % Full Newton-Raphson update
                        alpha = 1;
                    end
                    dU  = dU + alpha*ddU; % Updated solution
                    
                end
                
                % print info to screen
                fprintf(1, ['Load increment (',num2str(iincre),'/',num2str(nb_incre),')    NR iterations: ', num2str(i,'%3.0d'), '    relative residual: ', num2str(rel_res(i+1),'%1.2e'),'    time(s): ',num2str(toc(TimeIncr)) ,'\n']);
                
                Epl.bar = Epl.bar_trial;
                
            else
                break
            end
            
        end
        
    end
    
    
    
    if or(i==max_it_nr, dUmod==1) % did not reach convergence ---> try with displacement controlled
        
        dUmod = 1;
        
        dU    = dU_last; % reset dU
        
        %%% Reinitialization
        rel_res   	= NaN(1, max_it_nr+1);
        rel_res(1) 	= 1;
        
        bc_nod_ind = unique( MESH.SEGMENTS(:,MESH.segment_markers==6) );
        bc_ind     = [ ndim*(bc_nod_ind(:)-1)+1; ndim*(bc_nod_ind(:)-1)+2 ] ;
        bc_val     = dU(bc_ind);
        BC.ind     = [BC.ep_ind; bc_ind];
        BC.val     = [BC.ep_val; bc_val(:)];
        BC_NR       = BC;        % for Newton-Raphson correction calculation
        BC_NR.val   = 0*BC_NR.val; % for Newton-Raphson correction calculation
        
        for ii = 1:max_it_nr
            
            TimeIncr = tic;
            
            dEtrial      = strain_tot_calc(MESH, InvJ_all, dU, ieuv, nelblo);
            dEtrial.xy   = 2*dEtrial.xy; % gamma convention
            
            % Trial strain
            Etrial.xx = Etot.xx + dEtrial.xx;
            Etrial.yy = Etot.yy + dEtrial.yy;
            Etrial.zz = Etot.zz + dEtrial.zz;
            Etrial.xy = Etot.xy + dEtrial.xy;
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
            
            % Trial J2
            J2 = (Sdev.xx.^2 + Sdev.yy.^2 + Sdev.zz.^2)/2 + Sdev.xy.^2;
            
            % Check for plasticity consistency %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % cohe   = bsxfun(@plus, PARAM.ep_Co(MESH.elem_markers), PARAM.ep_H*Epl.bar);
            cohe   = bsxfun(@plus, PARAM.ep_Co(MESH.elem_markers) + PARAM.noise, PARAM.ep_H*Epl.bar);
            F      = sqrt(J2) + PARAM.ep_eta*P - PARAM.ep_xi*cohe;
            indx_p = F>0;
            
            d_gamma = zeros(nip,nel);
            
            if find(indx_p)
                
                % Return mapping with update of deviatoric stress tensor and P
%                 [d_g, Sdev, P, Epl, indx_apex] = Return_mapping_DruckerPrager(MESH, F, PARAM, indx_p, J2, Sdev, P, Epl);
                [d_g, Sdev, P, Epl, indx_apex] = Return_mapping_DruckerPrager_cohesion_noise(MESH, F, PARAM, indx_p, J2, Sdev, P, Epl);
                d_gamma(indx_p) = d_g;
                
                % Consistent tangent modulus (using non-corrected strain tensor!)
                D_el_ip = Tangent_Elastoplastic_DruckerPrager(MESH, indx_p, indx_apex, d_gamma, Etrial, PARAM, nip);
                
%                 % Update stiffness matrix
%                 if isequal(PARAM.ep_phi,PARAM.ep_xsi) % associative plasticity rule
                    elast = 0;
                    [A_all, ~, ~, ~] = matrix_entries_EP(MESH, D_el_ip, PARAM, SOLVER, elast);
%                 else
%                     Ael=elasticity_matrix(MESH.NODES,MESH.ELEMS,D_el_ip);
%                     A_all = reshape(Ael,size(Ael,1)*size(Ael,2),size(Ael,3));
%                 end
                
                
                % Compute residual
                Fint         = internal_force(MESH, InvJ_all, DetJ_all, Sdev, P + P_litho, SOLVER);
                res          =  - Fint;
                %                     opts_sparse.symmetric = 1;
                %                     opts_sparse.n_node_dof = 2;
                %                     A = sparse_create(MESH.ELEMS, A_all, opts_sparse);
                %                     opts.nthreads = [];
                %                     opts.symmetric  = 1;
                %                     Ac = sparse_convert(A, opts);
                %                     res = spmv(Ac, dU(:)) - rhs;
                res(BC_NR.ind)  = 0;
                rel_res(ii+1) = norm(res)/norm(Fext);
                
                if rel_res(ii+1)<rel_res_tol
                    % print info to screen
                    fprintf(1, ['Disp increment (',num2str(iincre),'/',num2str(nb_incre),')    NR iterations: ', num2str(ii,'%3.0d'), '    relative residual: ', num2str(rel_res(ii+1),'%1.2e'),'    time(s): ',num2str(toc(TimeIncr)) ,'\n\n']);
                    break
                else
                    % Calculate NR increment
                    ddU = mechanical_solver(MESH, A_all, BC_NR, res, SOLVER); % res as rhs
                    if isequal(SOLVER.ep_NR,'NR_line_search') % Line search
                        f = @(s) norm( res_calc_DP(MESH, InvJ_all, DetJ_all, Etot, Epl, PARAM, ieuv, BC_NR, indx_p, dU+s*ddU, Fext, P_litho, SOLVER) );
                        alpha = poly_search(f,norm(res));
                    elseif isequal(SOLVER.ep_NR,'NR_full_update') % Full Newton-Raphson update
                        alpha = 1;
                    end
                    dU  = dU + alpha*ddU; % Updated solution
                    
                end
                
                % print info to screen
                fprintf(1, ['Disp increment (',num2str(iincre),'/',num2str(nb_incre),')    NR iterations: ', num2str(ii,'%3.0d'), '    relative residual: ', num2str(rel_res(ii+1),'%1.2e'),'    time(s): ',num2str(toc(TimeIncr)) ,'\n']);
                
            else
                break
            end
            
        end
        
    end
    
    
    
    % Update converged solution
    Etot    = Etrial;
    U       = U + dU;
    dU_last = dU;
    Epl.bar = Epl.bar_trial;
    
    
    
    
    %%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if or(mod(iincre, POST.Hz_action)==0, iincre==1)
    if mod(iincre, POST.Hz_action)==0
             
        data_struct = [];
        xminp = POST.xminp;
        xmaxp = POST.xmaxp;
        yminp = POST.yminp;
        ymaxp = POST.ymaxp;
        
        if POST.action(4,1)==1
            figure(4); clf;
            if nnodel==3
                toplot = repmat(F,3,1);
            else
                toplot = repmat(mean(F),3,1);
            end
            trisurf(reshape(1:3*nel,3, nel)', MESH.NODES(1,MESH.ELEMS(1:3,:)), MESH.NODES(2,MESH.ELEMS(1:3,:)), zeros(size(MESH.NODES(1,MESH.ELEMS(1:3,:)))), toplot);
            title('Yield function F'); view(2); colorbar; axis image; axis off
            %             caxis([-1e8 1e9]);
            xlim([xminp xmaxp])
            ylim([yminp ymaxp])
            shading interp
            drawnow
            if POST.action(4,2)==1
                set(gcf,'PaperPositionMode','auto');
                print('-dpng',POST.fig_quality,[POST.folder_name,'/yield_function_i_',num2str(iincre,'%.4d')]);
            end
        end
        if POST.action(4,3)==1
            data_struct.yield_function = F;
        end
        
        if POST.action(5,1)==1
            figure(5); clf;
            if nnodel==3
                toplot = repmat(P+P_litho,3,1);
            else
                toplot = repmat(mean(P+P_litho),3,1);
            end
            trisurf(reshape(1:3*nel,3, nel)', MESH.NODES(1,MESH.ELEMS(1:3,:)), MESH.NODES(2,MESH.ELEMS(1:3,:)), zeros(size(MESH.NODES(1,MESH.ELEMS(1:3,:)))), toplot);
            title('Dynamic Pressure [Pa] (+ extension)'); view(2); colorbar; axis image; axis off
            caxis([min(toplot(:)) max(toplot(:))]);
            %             xmaxplot = max(abs(MESH.NODES(2,:)));
            %             xlim([-xmaxplot xmaxplot])
            xlim([xminp xmaxp])
            ylim([yminp ymaxp])
            shading interp
            drawnow
            if POST.action(5,2)==1
                set(gcf,'PaperPositionMode','auto');
                print('-dpng',POST.fig_quality,[POST.folder_name,'/solid_dyn_pressure_i_',num2str(iincre,'%.4d')]);
            end
        end
        if POST.action(5,3)==1
            data_struct.solid_dynamic_pressure = P+P_litho;
        end
        
       

        pSv=NaN(length(Sdev.xx),2);
        pSdmin=NaN(length(Sdev.xx),2);
        pSdmax=NaN(length(Sdev.xx),2);
%         S = zeros(2,2,2);

        for i = 1:length(Sdev.xx)
                mSxx = mean(Sdev.xx(:,i) + (P(:,i)));% + P_litho(:,i)));
                mSyy = mean(Sdev.yy(:,i) + (P(:,i)));% + P_litho(:,i)));
                mSxy = mean(Sdev.xy(:,i));
                S(:,:,1)=[mSxx,mSxy;mSxy,mSyy];
                
                [Vec,Val]=eig(S);
                I = sub2ind(size(Vec),[1 1],1:size(Vec,2));
                Vec = bsxfun(@times,Vec,sign(Vec(I)));
                Vec(:,2)   = -Vec(:,2) ;
                pSv(i,:)   = diag(Val);%.*sign(diag(prinstress));
                pSdmin(i,:)=Vec(:,1);
                pSdmax(i,:)=Vec(:,2);
        end
             
% figure(7),clf
% quiver(MESH.NODES(1,MESH.ELEMS(7,:)),MESH.NODES(2,MESH.ELEMS(7,:)), pSdmax(:,1)',pSdmax(:,2)','r')
% hold on
% quiver(MESH.NODES(1,MESH.ELEMS(7,:)),MESH.NODES(2,MESH.ELEMS(7,:)), pSdmin(:,1)',pSdmin(:,2)','b')

[X,Y] = meshgrid( linspace(min(MESH.NODES(1,:)),max(MESH.NODES(1,:)),800),  linspace(min(MESH.NODES(2,:)),max(MESH.NODES(2,:)),800) );
% Max_x = scatteredInterpolant(MESH.NODES(1,MESH.ELEMS(7,:))',MESH.NODES(2,MESH.ELEMS(7,:))',pSdmax(:,1)); 
% Max_y = scatteredInterpolant(MESH.NODES(1,MESH.ELEMS(7,:))',MESH.NODES(2,MESH.ELEMS(7,:))',pSdmax(:,2)); 
% pSdmax_X = Max_x(X,Y);
% pSdmax_Y = Max_y(X,Y);
% figure(7),clf,hold on
% % nb_s = 5;
% % pts_x = X(1:nb_s:end,1:nb_s:end); 
% % pts_y = Y(1:nb_s:end,1:nb_s:end);
%  intru_pts = MESH.NODES(:,MESH.node_markers==1);
intru_pts = MESH.NODES(:,MESH.node_markers==6);
pts_x = intru_pts(1,:);
pts_y = intru_pts(2,:);

rad = 1 + 0.0001 ;
theta = linspace(pi/2,-pi/2,30);
pts_x = rad*cos(theta) - 6.5;
pts_y = rad*sin(theta) ;

% rad = max(pts_y) + 0.0001 ;
% theta = linspace(pi/2,0,30);
% pts_x = (max(pts_x)/rad)*rad*cos(theta);
% pts_y = rad*sin(theta) ;

% h = streamline(X,Y,pSdmax_X,pSdmax_Y, pts_x, pts_y);
% set(h,'Color','r')
% drawnow
% Min_x = scatteredInterpolant(MESH.NODES(1,MESH.ELEMS(7,:))',MESH.NODES(2,MESH.ELEMS(7,:))',pSdmin(:,1)); 
% Min_y = scatteredInterpolant(MESH.NODES(1,MESH.ELEMS(7,:))',MESH.NODES(2,MESH.ELEMS(7,:))',pSdmin(:,2)); 
% pSdmin_X = Min_x(X,Y);
% pSdmin_Y = Min_y(X,Y);
% h = streamline(X,Y,pSdmin_X,pSdmin_Y, pts_x, pts_y);
% set(h,'Color','b')
% drawnow
% axis equal


% Maximum shear plane
pShearSd1 = 0*pSdmax;
pShearSd2 = 0*pSdmax;

fric_angle = 30; 

for i = 1:length(Sdev.xx)
    if MESH.NODES(2,MESH.ELEMS(7,i))<-1
        teta = -( 45 - fric_angle/2 ); 
        rotmat = [cosd(teta) -sind(teta); sind(teta) cosd(teta) ];
        pShearSd1(i,:) = pSdmin(i,:)*rotmat;
        teta = - (- 45 + fric_angle/2 ); 
        rotmat = [cosd(teta) -sind(teta); sind(teta) cosd(teta) ];
        pShearSd2(i,:) = pSdmin(i,:)*rotmat;
    else
        teta = -45 + fric_angle/2; 
        rotmat = [cosd(teta) -sind(teta); sind(teta) cosd(teta) ];
        pShearSd1(i,:) = pSdmin(i,:)*rotmat;
        teta = 45 - fric_angle/2; 
        rotmat = [cosd(teta) -sind(teta); sind(teta) cosd(teta) ];
        pShearSd2(i,:) = pSdmin(i,:)*rotmat;
    end
 
end
% figure(7),clf
% quiver(MESH.NODES(1,MESH.ELEMS(7,:)),MESH.NODES(2,MESH.ELEMS(7,:)), pSdmax(:,1)',pSdmax(:,2)','r')
% hold on
% quiver(MESH.NODES(1,MESH.ELEMS(7,:)),MESH.NODES(2,MESH.ELEMS(7,:)), pShearSd1(:,1)',pShearSd1(:,2)','b')
% quiver(MESH.NODES(1,MESH.ELEMS(7,:)),MESH.NODES(2,MESH.ELEMS(7,:)), pShearSd2(:,1)',pShearSd2(:,2)','g')
% axis equal

Max_x = scatteredInterpolant(MESH.NODES(1,MESH.ELEMS(7,:))',MESH.NODES(2,MESH.ELEMS(7,:))',pShearSd1(:,1)); 
Max_y = scatteredInterpolant(MESH.NODES(1,MESH.ELEMS(7,:))',MESH.NODES(2,MESH.ELEMS(7,:))',pShearSd1(:,2)); 
pShearSd1_X = Max_x(X,Y);
pShearSd1_Y = Max_y(X,Y);
Max_x = scatteredInterpolant(MESH.NODES(1,MESH.ELEMS(7,:))',MESH.NODES(2,MESH.ELEMS(7,:))',pShearSd2(:,1)); 
Max_y = scatteredInterpolant(MESH.NODES(1,MESH.ELEMS(7,:))',MESH.NODES(2,MESH.ELEMS(7,:))',pShearSd2(:,2)); 
pShearSd2_X = Max_x(X,Y);
pShearSd2_Y = Max_y(X,Y);
% figure(77),clf,hold on
% %  nb_s = 5;
% %  pts_x = X(1:nb_s:end,1:nb_s:end); 
% %  pts_y = Y(1:nb_s:end,1:nb_s:end);
% % % intru_pts = MESH.NODES(:,MESH.node_markers==6);
% % % pts_x = intru_pts(1,:);
% % % pts_y = intru_pts(2,:);
% h = streamline(X,Y,pShearSd1_X,pShearSd1_Y, pts_x, pts_y,[0.1 100]);
% set(h,'Color','g')
% hold on
% h = streamline(X,Y,pShearSd2_X,pShearSd2_Y, pts_x, pts_y,[0.1 100]);
% set(h,'Color','g')
% drawnow
% axis equal


% figure(8),clf
% toplot = repmat(squeeze(pSv(:,1)),3,1);
% trisurf(reshape(1:3*nel,3, nel)', MESH.NODES(1,MESH.ELEMS(1:3,:)), MESH.NODES(2,MESH.ELEMS(1:3,:)), zeros(size(MESH.NODES(1,MESH.ELEMS(1:3,:)))), toplot);
% view(2); colorbar
% shading interp
% 
% figure(9),clf
% toplot = repmat(squeeze(pSv(:,2)),3,1);
% trisurf(reshape(1:3*nel,3, nel)', MESH.NODES(1,MESH.ELEMS(1:3,:)), MESH.NODES(2,MESH.ELEMS(1:3,:)), zeros(size(MESH.NODES(1,MESH.ELEMS(1:3,:)))), toplot);
% view(2); colorbar
% shading interp 



        if POST.action(6,1)==1
            figure(6); clf;
            subplot(221)
            if nnodel==3
                toplot = repmat(Sdev.xx,3,1);
            else
                toplot = repmat(mean(Sdev.xx),3,1);
            end
            trisurf(reshape(1:3*nel,3, nel)', MESH.NODES(1,MESH.ELEMS(1:3,:)), MESH.NODES(2,MESH.ELEMS(1:3,:)), zeros(size(MESH.NODES(1,MESH.ELEMS(1:3,:)))), toplot);
            title('Dev Stress xx [Pa]'); view(2); colorbar; axis image; axis off
            caxis([min(toplot(:)) max(toplot(:))]);
            xlim([xminp xmaxp])
            ylim([yminp ymaxp])
            shading interp
            subplot(2,2,2)
            if nnodel==3
                toplot = repmat(Sdev.yy,3,1);
            else
                toplot = repmat(mean(Sdev.yy),3,1);
            end
            trisurf(reshape(1:3*nel,3, nel)', MESH.NODES(1,MESH.ELEMS(1:3,:)), MESH.NODES(2,MESH.ELEMS(1:3,:)), zeros(size(MESH.NODES(1,MESH.ELEMS(1:3,:)))), toplot);
            title('Dev Stress yy [Pa]'); view(2); colorbar; axis image; axis off
            caxis([min(toplot(:)) max(toplot(:))]);
            shading interp
            subplot(2,2,[3,4])
            if nnodel==3
                toplot = repmat(Sdev.xy,3,1);
            else
                toplot = repmat(mean(Sdev.xy),3,1);
            end
            trisurf(reshape(1:3*nel,3, nel)', MESH.NODES(1,MESH.ELEMS(1:3,:)), MESH.NODES(2,MESH.ELEMS(1:3,:)), zeros(size(MESH.NODES(1,MESH.ELEMS(1:3,:)))), toplot);
            title('Dev Stress xy [Pa] (mean value per element)'); view(2); colorbar; axis image; axis off
            caxis([min(toplot(:)) max(toplot(:))]);
            xlim([xminp xmaxp])
            ylim([yminp ymaxp])
            shading interp
            if POST.action(6,2)==1
                set(gcf,'PaperPositionMode','auto');
                print('-dpng',POST.fig_quality,[POST.folder_name,'/dev_stress_i_',num2str(iincre,'%.4d')]);
            end
        end
        if POST.action(6,3)==1
            data_struct.dev_stress = Sdev;
        end
        
        if POST.action(7,1)==1
            figure(7),clf
            E2 = sqrt( ( (Etot.xx- Etrial.ii/3).^2 + (Etot.yy- Etrial.ii/3).^2 + (Etot.zz- Etrial.ii/3).^2 )/2 + (Etrial.xy/2).^2 );
            if nnodel==3
                toplot = repmat(E2,3,1);
            else
                toplot = repmat(mean(E2),3,1);
            end
            trisurf(reshape(1:3*nel,3, nel)', MESH.NODES(1,MESH.ELEMS(1:3,:)), MESH.NODES(2,MESH.ELEMS(1:3,:)), zeros(size(MESH.NODES(1,MESH.ELEMS(1:3,:)))), toplot);
            title('Dev Strain 2nd Inv (mean value per element)'); view(2); colorbar; axis image; axis off
            xlim([xminp xmaxp])
            ylim([yminp ymaxp])
            caxis([min(toplot(:)) max(toplot(:))]);
            shading interp
            drawnow
            if POST.action(7,2)==1
                set(gcf,'PaperPositionMode','auto');
                print('-dpng',POST.fig_quality,[POST.folder_name,'/dev_strain_2nd_i_',num2str(iincre,'%.4d')]);
            end
        end
        if POST.action(7,3)==1
            data_struct.tot_strain = Etot;
        end
        
        if POST.action(8,1)==1
            figure(8),clf
            if nnodel==3
                toplot = repmat(Etrial.ii,3,1);
            else
                toplot = repmat(mean(Etrial.ii),3,1);
            end
            trisurf(reshape(1:3*nel,3, nel)', MESH.NODES(1,MESH.ELEMS(1:3,:)), MESH.NODES(2,MESH.ELEMS(1:3,:)), zeros(size(MESH.NODES(1,MESH.ELEMS(1:3,:)))), toplot);
            title('Vol Strain (mean value per element)'); view(2); colorbar; axis image; axis off
            xlim([xminp xmaxp])
            ylim([yminp ymaxp])
            caxis([min(toplot(:)) max(toplot(:))]);
            shading interp
            drawnow
            if POST.action(8,2)==1
                set(gcf,'PaperPositionMode','auto');
                print('-dpng',POST.fig_quality,[POST.folder_name,'/vol_strain_i_',num2str(iincre,'%.4d')]);
            end
        end
        if POST.action(8,3)==1
            data_struct.vol_strain = Etot;
        end
        
        if POST.action(9,1)==1
            figure(9),clf
            if nnodel==3
                toplot = repmat(Epl.bar,3,1);
            else
                toplot = repmat(mean(Epl.bar),3,1);
            end
            
            toplot(toplot<1e-10) = nan;
            
            trisurf(reshape(1:3*nel,3, nel)', MESH.NODES(1,MESH.ELEMS(1:3,:)), MESH.NODES(2,MESH.ELEMS(1:3,:)), zeros(size(MESH.NODES(1,MESH.ELEMS(1:3,:)))), toplot);
            title('Acc. Plastic Strain (mean value per element)'); view(2); colorbar;
            
            if isnan(toplot<1e-10) 
            caxis([min(toplot(:)) max(toplot(:))]);
            end
            %             xmaxplot = max(abs(MESH.NODES(2,:)));
            %             xlim([-xmaxplot xmaxplot])
            %             xmaxplot = 10000;
            %             xlim([0  xmaxplot])
            %             ylim([-xmaxplot 0])

             
            hold on
            h = streamline(X,Y,pShearSd1_X,pShearSd1_Y, pts_x, pts_y,[.5 170]);
            set(h,'Color','k')
            hold on
            h = streamline(X,Y,pShearSd2_X,pShearSd2_Y, pts_x, pts_y,[.5 170]);
            set(h,'Color','k')
            
            
            xlim([xminp xmaxp])
            ylim([yminp ymaxp])
            
            shading interp
            box on
            grid on
            
            rad = 1 ;
            theta = linspace(pi/2,0,30);
            pts_x = rad*cos(theta);
            pts_y = rad*sin(theta) ;
            hold on
            plot(pts_x, pts_y,'-k')
            drawnow
            
            if POST.action(9,2)==1
                set(gcf,'PaperPositionMode','auto');
                print('-dpng',POST.fig_quality,[POST.folder_name,'/acc_plast_strain_i_',num2str(iincre,'%.4d')]);
            end
        end
        if POST.action(9,3)==1
            data_struct.plast_strain = Epl.bar;
        end
        
        if POST.action(10,1)==1
            figure(10),clf
            subplot(221)
            toplot = U(1:2:end-1) ;
            toplot = toplot(MESH.ELEMS(1:3,:));
            trisurf(reshape(1:3*nel,3, nel)', MESH.NODES(1,MESH.ELEMS(1:3,:)), MESH.NODES(2,MESH.ELEMS(1:3,:)), zeros(size(MESH.NODES(1,MESH.ELEMS(1:3,:)))), toplot);
            title('U xx'); view(2); colorbar; axis image; shading interp 
            xlim([xminp xmaxp])
            ylim([yminp ymaxp])
            caxis([min(toplot(:)) max(toplot(:))]);
            shading interp
            subplot(2,2,2)
            toplot = U(2:2:end);
            toplot = toplot(MESH.ELEMS(1:3,:));
            trisurf(reshape(1:3*nel,3, nel)', MESH.NODES(1,MESH.ELEMS(1:3,:)), MESH.NODES(2,MESH.ELEMS(1:3,:)), zeros(size(MESH.NODES(1,MESH.ELEMS(1:3,:)))), toplot);
            title('U yy'); view(2); colorbar; axis image; axis off
            xlim([xminp xmaxp])
            ylim([yminp ymaxp])
            caxis([min(toplot(:)) max(toplot(:))]);
            shading interp
            subplot(2,2,[3,4]), hold on
            surf_nod = unique(MESH.SEGMENTS(:, find(MESH.segment_markers==3)));
            toplot = U(2:2:end);
            plot(MESH.NODES(1,surf_nod),toplot(surf_nod),'rd')
            plot([min(MESH.NODES(1,:)) max(MESH.NODES(1,:))] , [0 0],'-k')
            xlim([xminp xmaxp])
            %             ylim([yminp ymaxp])
            title('Surface vertical displacement [m]')
            drawnow
            if POST.action(9,2)==1
                set(gcf,'PaperPositionMode','auto');
                print('-dpng',POST.fig_quality,[POST.folder_name,'/solid_disp_i_',num2str(iincre,'%.4d')]);
            end
        end
        if POST.action(10,3)==1
            data_struct.solid_disp = U;
        end
        
        if isempty(data_struct)==0
            save([POST.folder_name,'/data_i_',num2str(iincre,'%.4d')],'data_struct')
        end
        
    end
    
    
    
%      
%     figure(11),clf
%     
%     subplot(211)
%     if nnodel==3
%         toplot = repmat(Epl.bar,3,1);
%     else
%         toplot = repmat(mean(Epl.bar),3,1);
%     end
%     trisurf(reshape(1:3*nel,3, nel)', MESH.NODES(1,MESH.ELEMS(1:3,:)), MESH.NODES(2,MESH.ELEMS(1:3,:)), zeros(size(MESH.NODES(1,MESH.ELEMS(1:3,:)))), toplot);
%     title('Acc. Plastic Strain (mean value per element)'); view(2); colorbar; axis image; axis off
%     caxis([min(toplot(:)) max(toplot(:))]);
%     if max(toplot(:))>1
%         caxis([min(toplot(:)) 1]);
%     end
%     xlim([xminp xmaxp])
%     ylim([yminp ymaxp])
%     shading interp
%     drawnow
%     
%     subplot(212)
%     if nnodel==3
%         toplot = repmat(F,3,1);
%     else
%         toplot = repmat(mean(F),3,1);
%     end
%     trisurf(reshape(1:3*nel,3, nel)', MESH.NODES(1,MESH.ELEMS(1:3,:)), MESH.NODES(2,MESH.ELEMS(1:3,:)), zeros(size(MESH.NODES(1,MESH.ELEMS(1:3,:)))), toplot);
%     title('Yield function F'); view(2); colorbar; axis image; axis off
%     xlim([xminp xmaxp])
%     ylim([yminp ymaxp])
%     shading interp
%     drawnow
%     
%     if POST.action(1,2)==1
%     set(gcf,'PaperPositionMode','auto');
%     print('-dpng','-r120',[POST.folder_name,'/anim_i_',num2str(iincre,'%.4d')]);
%     end
    
    
    
    
    
    
    
%     
%     %%% TESTING ELASTICITY SURF DEFORMATION WITH McTigue 1988
%     
%     xi = linspace(-10000, 10000, 500);
%     yi = 0*xi;
%     zi = 0*xi;
%     
%     addpath(genpath('/Users/albanes/Desktop/DIPS/Modelling/USGS TM13'));
%     P_G = (BC.v_press / SOLVER.ep_nb_incre) * iincre / PARAM.ep_G(1);
%     nu = PARAM.ep_nu;
%     a = 500;
%     z0 = 2000;
%     [u_a v_a w_a dwdx_a dwdy_a] = mctigue2D(0,0,z0,P_G,a,nu,xi,yi);   % 2D analytical model
%     
%     fix = scatteredInterpolant(MESH.NODES(1,:)',MESH.NODES(2,:)',MESH.NODES(3,:)', U(1:3:end-2) );
%     fiy = scatteredInterpolant(MESH.NODES(1,:)',MESH.NODES(2,:)',MESH.NODES(3,:)', U(2:3:end-1) );
%     fiz = scatteredInterpolant(MESH.NODES(1,:)',MESH.NODES(2,:)',MESH.NODES(3,:)', U(3:3:end) );
%     ux_fem = fix(xi, yi,zi);
%     uy_fem = fiy(xi, yi,zi);
%     uz_fem = fiz(xi, yi,zi);
%     
%     figure(111),clf
%     subplot(131), hold on
%     plot(xi, ux_fem, 'g', xi, u_a,'sr'), ylabel('horizontal EW displacement')
%     subplot(132), hold on
%     plot(xi, uy_fem, 'g', xi, v_a,'sr'), ylabel('horizontal NS displacement')
%     subplot(133), hold on
%     plot(xi, uz_fem, 'g', xi, w_a,'sr'), ylabel('vertical displacement')
%     
%     drawnow
%         
    
    
    
end






