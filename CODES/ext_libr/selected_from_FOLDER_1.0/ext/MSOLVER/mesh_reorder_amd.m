function [MESH, perm, iperm] = mesh_reorder_amd(MESH)

%% NODES
spopts.symmetric = 0;
spopts.ndof      = 1;
A = sparse_create(MESH.ELEMS, 1, spopts);
perm = uint32(amd(A));
iperm = zeros(size(perm), 'uint32');
iperm(perm)  = 1:length(perm);

MESH.ELEMS = reshape(iperm(MESH.ELEMS), size(MESH.ELEMS));
MESH.NODES = MESH.NODES(:,perm);

if isfield(MESH, 'node_markers')
    MESH.node_markers = MESH.node_markers(:,perm);
end
if isfield(MESH, 'EDGES')
    MESH.EDGES = iperm(MESH.EDGES);
end
if isfield(MESH, 'SEGMENTS')
    MESH.SEGMENTS = iperm(MESH.SEGMENTS);
end

%% ELEMS
[~, permel]  = sort(min(MESH.ELEMS));
permel       = uint32(permel);
ipermel      = zeros(size(permel), 'uint32');
ipermel(permel) = 1:length(permel);

MESH.ELEMS  = MESH.ELEMS(:,permel);
if isfield(MESH, 'elem_markers')
    MESH.elem_markers = MESH.elem_markers(permel);
end
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
