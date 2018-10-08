function MESH = NormalEdges( MESH )

% Part of Open-GeoNabla, copyright GPLv3, 2018
% https://github.com/albansouche/Open-GeoNabla/
% Physics of Geological Processes (PGP) , The NJORD Centre, Dept of Geosciences, University of Oslo
% Author: Alban Souche


nnodel = size(MESH.ELEMS,1);
normal=[-( MESH.NODES(2, MESH.SEGMENTS(2,:)) - MESH.NODES(2, MESH.SEGMENTS(1,:)) ); ...
    MESH.NODES(1, MESH.SEGMENTS(2,:)) - MESH.NODES(1, MESH.SEGMENTS(1,:)) ];
midpoints = (MESH.NODES(:,MESH.SEGMENTS(1,:)) + MESH.NODES(:,MESH.SEGMENTS(2,:)) )./2;
if nnodel==7
outwardvector= midpoints - MESH.NODES(:, MESH.ELEMS(7,MESH.elem_segments)); 
elseif nnodel==3
outwardvector= midpoints - ( MESH.NODES(:, MESH.ELEMS(1,MESH.elem_segments))+MESH.NODES(:, MESH.ELEMS(2,MESH.elem_segments))+MESH.NODES(:, MESH.ELEMS(3,MESH.elem_segments)) )./3  ;     
end
%check whether normal is oriented in the outward direction
temp=dot(normal, outwardvector)<0;
normal(:,temp)=-normal(:,temp);
MESH.normal=normal./repmat(abs(normal(1, :) + 1i*normal(2, :)), 2, 1);


% % Original from Jan
% [ndim, ~]=size(MESH.NODES);
% nseg=size(MESH.SEGMENTS, 2);
% 
% outwardvector=zeros(ndim, nseg);
% midpoints=zeros(ndim, nseg);
% 
% normal=[-( MESH.NODES(2, MESH.SEGMENTS(2,:)) - MESH.NODES(2, MESH.SEGMENTS(1,:)) ); ...
%     MESH.NODES(1, MESH.SEGMENTS(2,:)) - MESH.NODES(1, MESH.SEGMENTS(1,:)) ];
% 
% midpoints(1, :)=(MESH.NODES(1, MESH.SEGMENTS(1, :)) + MESH.NODES(1, MESH.SEGMENTS(2, :)))/2;
% midpoints(2, :)=(MESH.NODES(2, MESH.SEGMENTS(1, :)) + MESH.NODES(2, MESH.SEGMENTS(2, :)))/2;
% 
% ind_needed=zeros(1, nseg);
% elems=MESH.ELEMS(1:3, :);
% for i=1:nseg
%     ind=( elems(:, MESH.elem_segments(i))==MESH.SEGMENTS(1, i))...
%         + (elems(:, MESH.elem_segments(i))==MESH.SEGMENTS(2, i) );
%     
%     ind_needed(i)= elems(~ind, MESH.elem_segments(i));
% end
% 
% outwardvector(1, :)= midpoints(1, :) - MESH.NODES(1, ind_needed);
% outwardvector(2, :)= midpoints(2, :) - MESH.NODES(2, ind_needed);
% 
% %check whether normal is oriented in the outward direction
% temp=dot(normal, outwardvector)<0;
% normal(:,temp)=-normal(:,temp);
% 
% MESH.normal=normal./repmat(abs(normal(1, :) + 1i*normal(2, :)), 2, 1);
