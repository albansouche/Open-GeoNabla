function P_litho = initial_litho_zlevel(MESH, PARAM, SOLVER, POST)

% Part of Open-GeoNabla, copyright GPLv3, 2018
% https://github.com/albansouche/Open-GeoNabla/
% Physics of Geological Processes (PGP) , The NJORD Centre, Dept of Geosciences, University of Oslo
% Author: Alban Souche

if ~isequal(PARAM.g(2),0)

    [nnodel,nel] = size(MESH.ELEMS);
    
    %% Integration and evaluation points
    neva      = SOLVER.neva;
    [ieuv, ~] = ip_triangle(neva);
    P_zlevel =  PARAM.g(2)*PARAM.ep_rho(1)*max(MESH.NODES(2,:));
    P_nod     =  P_zlevel * ones(size(MESH.NODES(2,:)));
    for i = 1:max(MESH.elem_markers)
        nod_i = MESH.ELEMS(:,MESH.elem_markers==i);
        ref_z = max(MESH.NODES(2,nod_i));
        max_P = max(P_nod);
        P_nod(nod_i) = max_P + PARAM.g(2)*PARAM.ep_rho(i)*(MESH.NODES(2,nod_i)-ref_z);
    end
    
    % Pressure Calculation at integration points
    nelblo       = min(nel, SOLVER.nelblo);
    nblo         = ceil(nel/nelblo);
    il           = 1;
    iu           = nelblo;
    % Derivative of shape functions
    [ N, ~]   = shp_deriv_triangle(ieuv, nnodel);
    % Allocation
    P_litho = zeros(neva,nel);
    for ib = 1:nblo
        % Extract block Pressure
        P_block = reshape(P_nod(MESH.ELEMS(:,il:iu)), nnodel, nelblo);
        for ieva=1:neva
            %Extract spatial derivative
            Ni = N{ieva};
            %Calculate Pressure at itegration point
            P_litho(ieva,il:iu) = Ni'*P_block;
        end
        il  = il+nelblo;
        if(ib==nblo-1)
            nelblo 	  = nel-iu;
        end
        iu  = iu+nelblo;
    end
    
    if POST.action(5,1)==1
        figure(5); clf;
        if nnodel==3
            toplot = repmat(P_litho,3,1);
        else
            toplot = repmat(mean(P_litho),3,1);
        end
        trisurf(reshape(1:3*nel,3, nel)', MESH.NODES(1,MESH.ELEMS(1:3,:)), MESH.NODES(2,MESH.ELEMS(1:3,:)), zeros(size(MESH.NODES(1,MESH.ELEMS(1:3,:)))), toplot);
        title('Lithostatic Pressure [Pa] (mean value per element)'); view(2); colorbar; axis image; axis off
        caxis([min(toplot(:)) max(toplot(:))]);
        xlim([POST.xminp POST.xmaxp])
        ylim([POST.yminp POST.ymaxp])
        shading interp
        drawnow
        if POST.action(5,2)==1
            set(gcf,'PaperPositionMode','auto');
            print('-dpng',POST.fig_quality,[POST.folder_name,'/solid_litho_pressure_i_',num2str(1,'%.4d')]); 
        end
    end
    
    
else
    
    P_litho        = 0;

end