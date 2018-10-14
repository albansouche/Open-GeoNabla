function MESH = Normal_calc_2D(MESH)

% Part of Open-GeoNabla, copyright GPLv3, 2018
% https://github.com/albansouche/Open-GeoNabla/
% Physics of Geological Processes (PGP) , The NJORD Centre, Dept of Geosciences, University of Oslo
% Author: Alban Souche


if size(MESH.ELEMS,1)==7
    MESH.SEGMENTS = [MESH.SEGMENTS; ones(1,length(MESH.SEGMENTS))];
end

MESH.elem_segments = zeros(1,length(MESH.SEGMENTS));

% Set important tsearch2 parameters
WS = [];
WS.NEIGHBORS = MESH.NEIGHBORS;  % element neighbor information
WS.xmin = min(MESH.NODES(1,:)); % domain extents
WS.xmax = max(MESH.NODES(1,:));
WS.ymin = min(MESH.NODES(2,:));
WS.ymax = max(MESH.NODES(2,:));


% Compute the normal to the facets %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
normal=[-( MESH.NODES(2, MESH.SEGMENTS(2,:)) - MESH.NODES(2, MESH.SEGMENTS(1,:)) ); ...
    MESH.NODES(1, MESH.SEGMENTS(2,:)) - MESH.NODES(1, MESH.SEGMENTS(1,:)) ];
MESH.normal=normal./repmat(abs(normal(1, :) + 1i*normal(2, :)), 2, 1);

% Run tsearch2 on 2 CPUs
opts.nthreads = 2;
markers       = (MESH.NODES(:,MESH.SEGMENTS(1,:)) + MESH.NODES(:,MESH.SEGMENTS(2,:)) )./2;
% check on left and right of the segments to find the corresponding element
elem_plus_side  = tsearch2(MESH.NODES, MESH.ELEMS(1:3, :), markers + (1e-6.*MESH.normal), WS, [], opts);
elem_minus_side = tsearch2(MESH.NODES, MESH.ELEMS(1:3, :), markers - (1e-6.*MESH.normal), WS, [], opts);
MESH.elem_segments = max(elem_plus_side,elem_minus_side);    % the normal vector should be pointing outward along boundary facets
MESH.normal(:,elem_plus_side>0) = -MESH.normal(:,elem_plus_side>0); % the normal vector should be pointing outward along boundary facets

[~, UVm] = einterp(MESH, ones(1,length(MESH.NODES)), markers, MESH.elem_segments, opts);
vec_nod2find = (UVm<1e-6);
vec_nod4     = false(size(vec_nod2find)); 
vec_nod5     = [true(1,length(vec_nod2find)); false(1,length(vec_nod2find))];  
vec_nod6     = [false(1,length(vec_nod2find)); true(1,length(vec_nod2find))];  

if size(MESH.ELEMS,1)>3
    indx = sum(vec_nod2find==vec_nod4)==2;
    MESH.SEGMENTS(3,indx) = MESH.ELEMS(4,MESH.elem_segments(indx));
    indx = sum(vec_nod2find==vec_nod5)==2;
    MESH.SEGMENTS(3,indx) = MESH.ELEMS(5,MESH.elem_segments(indx));
    indx = sum(vec_nod2find==vec_nod6)==2;
    MESH.SEGMENTS(3,indx) = MESH.ELEMS(6,MESH.elem_segments(indx));
end

