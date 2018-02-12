function P_litho = initial_litho(MESH, PARAM, SOLVER, POST)

if ~isequal(PARAM.g(3),0)

    [nnodel,nel] = size(MESH.ELEMS);
    
    %% Integration and evaluation points
    neva      = SOLVER.neva;
    [ieuv, ~] = ip_tetra(neva);
    P_nod     = zeros(size(MESH.NODES(3,:)));
    for i = 1:max(MESH.elem_markers)
        nod_i = MESH.ELEMS(:,MESH.elem_markers==i);
        ref_z = max(MESH.NODES(3,nod_i));
        max_P = max(P_nod);
        P_nod(nod_i) = max_P + PARAM.g(3)*PARAM.ep_rho(i)*(MESH.NODES(3,nod_i)-ref_z);
    end
    
    % Pressure Calculation at integration points
    nelblo       = min(nel, SOLVER.nelblo);
    nblo         = ceil(nel/nelblo);
    il           = 1;
    iu           = nelblo;
    % Derivative of shape functions
    [ N, ~]   = shp_deriv_tetra(ieuv, nnodel);
    
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
    
if POST.action(3)==1
    to_data_n = zeros(1,size(MESH.NODES,2));
    elems_tmp = MESH.ELEMS(1:4,:);
    if size(MESH.ELEMS,1)==10
        to_data_n(elems_tmp(:)) = P_litho(:);
    elseif size(MESH.ELEMS,1)==4
        to_data_n(elems_tmp(:)) = repmat(mean(P_litho,1),4,1);
    end
    paraview_write('tet4' ,MESH.NODES, MESH.ELEMS, MESH.elem_markers, to_data_n, [POST.folder_name,'/P_litho'], 0)
end
    
    
else
    
    P_litho        = 0;

end