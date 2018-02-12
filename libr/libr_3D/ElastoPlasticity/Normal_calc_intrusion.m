function MESH2 = Normal_calc_intrusion(MESH, MESHARG, SOLVER)

% By ALBAN SOUCHE, Physics of Geological Prcesses, Njord Center, Departement of Geosciences, University of Oslo, January 2017

% Save MESH in MESH2
MESH2         = MESH;

% Only compute the normals of the intrusion facets
face_id       = MESHARG.face_id(7);
indx_face_out = ~(MESH.face_markers==face_id);
MESH.FACES(:,indx_face_out)      = []; % delete facets except of intrusion facets

% Compute the normal to the facets %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U = [ MESH.NODES(1, MESH.FACES(2,:)) - MESH.NODES(1, MESH.FACES(1,:)) ; ...
      MESH.NODES(2, MESH.FACES(2,:)) - MESH.NODES(2, MESH.FACES(1,:)) ; ...
      MESH.NODES(3, MESH.FACES(2,:)) - MESH.NODES(3, MESH.FACES(1,:)) ];
V = [ MESH.NODES(1, MESH.FACES(3,:)) - MESH.NODES(1, MESH.FACES(1,:)) ; ...
      MESH.NODES(2, MESH.FACES(3,:)) - MESH.NODES(2, MESH.FACES(1,:)) ; ...
      MESH.NODES(3, MESH.FACES(3,:)) - MESH.NODES(3, MESH.FACES(1,:)) ];
normal =[ U(2,:).*V(3,:) - U(3,:).*V(2,:) ; ...
          U(3,:).*V(1,:) - U(1,:).*V(3,:) ; ...
          U(1,:).*V(2,:) - U(2,:).*V(1,:) ];
normal = bsxfun(@rdivide,  normal , sqrt(normal(1,:).^2 + normal(2,:).^2 + normal(3,:).^2) ) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% Find the element of the facets %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MESH.elem_faces = zeros(1,length(MESH.FACES));

% Set important tsearch2 parameters
WS = [];
WS.NEIGHBORS = uint32(MESH.NEIGHBORS);  % element neighbor information
WS.xmin = min(MESH.NODES(1,:)); % domain extents
WS.xmax = max(MESH.NODES(1,:));
WS.ymin = min(MESH.NODES(2,:));
WS.ymax = max(MESH.NODES(2,:));
WS.zmin = min(MESH.NODES(3,:));
WS.zmax = max(MESH.NODES(3,:));

opts.nthreads = SOLVER.nb_cpu;
mid_pts       = (MESH.NODES(:,MESH.FACES(1,:)) + MESH.NODES(:,MESH.FACES(2,:)) + MESH.NODES(:,MESH.FACES(3,:)) )./3;
% check on left and right of the segments to find the corresponding element
elem_plus_side  = tsearch2(MESH.NODES, MESH.ELEMS(1:4, :), mid_pts + (1e-6.*normal), WS, [], opts);

if find(elem_plus_side==0)
    elem_minus_side = tsearch2(MESH.NODES, MESH.ELEMS(1:4, :), mid_pts - (1e-6.*normal), WS, [], opts);
else
    elem_minus_side = 0*elem_plus_side;
end

if find(MESH.elem_markers(elem_plus_side)==6)
    elem_plus_side(MESH.elem_markers(elem_plus_side)==6)   = 0; % if elements is part of the intrusion (phase=6), then 0. We only want element number from the solid mesh.
    elem_minus_side(MESH.elem_markers(elem_minus_side)==6) = 0; % if elements is part of the intrusion (phase=6), then 0. We only want element number from the solid mesh.
end

MESH.elem_faces = max(elem_plus_side,elem_minus_side);
normal(:,elem_plus_side>0) = -normal(:,elem_plus_side>0); % the normal vector should be pointing outward along boundary facets (in respect to the solid mesh)

MESH2.normal = zeros(size(MESH2.NODES,1), size(MESH2.FACES,2));
MESH2.normal(:,MESH.face_markers==face_id) = normal;
