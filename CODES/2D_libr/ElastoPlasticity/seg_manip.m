function MESH = seg_manip(MESH,type_el, SOLVER)

% Part of Open-GeoNabla, copyright GPLv3, 2018
% https://github.com/albansouche/Open-GeoNabla/
% Physics of Geological Processes (PGP) , The NJORD Centre, Dept of Geosciences, University of Oslo
% Author: Alban Souche


if isequal(type_el,'tri7')
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

opts.nthreads = SOLVER.nb_cpu;
mid_pts       = (MESH.NODES(:,MESH.SEGMENTS(1,:)) + MESH.NODES(:,MESH.SEGMENTS(2,:)) )./2;
MESH.elem_segments = tsearch2(MESH.NODES, MESH.ELEMS(1:3, :), mid_pts, WS, [], opts);

% in the case tsearch2 return some 0 value... check on left and right of
% the segments to find the corresponding element
indx0 = find(MESH.elem_segments==0);
if indx0
    mid_pts_right  = [ mid_pts(1,indx0) - (MESH.NODES(2,MESH.SEGMENTS(2,indx0)) - MESH.NODES(2,MESH.SEGMENTS(1,indx0)) )./1e3; ...
        mid_pts(2,indx0) + (MESH.NODES(1,MESH.SEGMENTS(2,indx0)) - MESH.NODES(1,MESH.SEGMENTS(1,indx0)) )./1e3 ] ;
    mid_pts_left = [ mid_pts(1,indx0) + (MESH.NODES(2,MESH.SEGMENTS(2,indx0)) - MESH.NODES(2,MESH.SEGMENTS(1,indx0)) )./1e3; ...
        mid_pts(2,indx0) - (MESH.NODES(1,MESH.SEGMENTS(2,indx0)) - MESH.NODES(1,MESH.SEGMENTS(1,indx0)) )./1e3 ] ;
    elem_left  = tsearch2(MESH.NODES, MESH.ELEMS(1:3, :), mid_pts_left, WS, [], opts);
    elem_right = tsearch2(MESH.NODES, MESH.ELEMS(1:3, :), mid_pts_right, WS, [], opts);
    MESH.elem_segments(indx0) = max(elem_left,elem_right);
end

if isequal(type_el,'tri7')
    [~, UVm] = einterp(MESH, ones(1,length(MESH.NODES)), mid_pts, MESH.elem_segments, opts);
    vec_nod2find = (UVm<1e-6);
    vec_nod4     = false(size(vec_nod2find));
    vec_nod5     = [true(1,length(vec_nod2find)); false(1,length(vec_nod2find))];
    vec_nod6     = [false(1,length(vec_nod2find)); true(1,length(vec_nod2find))];
    
    indx = sum(vec_nod2find==vec_nod4)==2;
    MESH.SEGMENTS(3,indx) = MESH.ELEMS(4,MESH.elem_segments(indx));
    indx = sum(vec_nod2find==vec_nod5)==2;
    MESH.SEGMENTS(3,indx) = MESH.ELEMS(5,MESH.elem_segments(indx));
    indx = sum(vec_nod2find==vec_nod6)==2;
    MESH.SEGMENTS(3,indx) = MESH.ELEMS(6,MESH.elem_segments(indx));
end


% % % Old from Jan's routine
% % for i=1:length(MESH.SEGMENTS)
% %     tmp = ismember(MESH.ELEMS,MESH.SEGMENTS(1:2,i));
% %     ele_seg = find(sum(tmp)==2);
% %     
% %     if isequal(type_el,'tri7');
% %         nod_loc = tmp(:,ele_seg)';
% %         if isequal(nod_loc(1:3),[1 1 0])
% %             MESH.SEGMENTS(3,i)=MESH.ELEMS(6,ele_seg);
% %         elseif isequal(nod_loc(1:3),[1 0 1])
% %             MESH.SEGMENTS(3,i)=MESH.ELEMS(5,ele_seg);
% %         elseif isequal(nod_loc(1:3),[0 1 1])
% %             MESH.SEGMENTS(3,i)=MESH.ELEMS(4,ele_seg);
% %         end
% %     end
% %     MESH.elem_segments(i) = ele_seg;
% % end



