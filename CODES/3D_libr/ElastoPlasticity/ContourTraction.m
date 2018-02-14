function [F, cont_tx, cont_ty, cont_tz] = ContourTraction(MESH, SIGMA, seg_id)

% Part of Open-GeoNabla, copyright GPLv3, 2018
% https://github.com/albansouche/Open-GeoNabla/
% Physics of Geological Processes (PGP) , The NJORD Centre, Dept of Geosciences, University of Oslo
% Author: Jan Cornet

% ContourTraction computes the contour integral of tractions 

NODES       = MESH.NODES;
[ndim,nnod] = size(NODES);

% Traction only along specific segment number %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
face_id = seg_id;
indx_face_out = ~(MESH.face_markers==face_id);
MESH.face_markers(indx_face_out) = [];
MESH.FACES(:,indx_face_out)      = [];
MESH.normal(:,indx_face_out)     = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SEGFA=MESH.FACES;
[nbound,nn]=size(SEGFA);
ipw  = 0.5;
N    = ones(3,1)/3;
dNdU = [-1 , -1 ; 1 , 0 ; 0, 1];

sigma=zeros(nbound, nn, ndim*(ndim+1)/2);
sigma(:, :, 1)=reshape(SIGMA(1, SEGFA), size(SEGFA));
sigma(:, :, 2)=reshape(SIGMA(2, SEGFA), size(SEGFA));
sigma(:, :, 3)=reshape(SIGMA(3, SEGFA), size(SEGFA));

TT=zeros(nbound,nn,ndim);
sigma(:, :, 4)=reshape(SIGMA(4, SEGFA), size(SEGFA));
sigma(:, :, 5)=reshape(SIGMA(5, SEGFA), size(SEGFA));
sigma(:, :, 6)=reshape(SIGMA(6, SEGFA), size(SEGFA));

TT(:,:,1)=sigma(:, :, 1).*repmat(MESH.normal(1, :), nbound, 1) + sigma(:, :, 4).*repmat(MESH.normal(2, :), nbound, 1) + sigma(:, :, 5).*repmat(MESH.normal(3, :), nbound, 1);
TT(:,:,2)=sigma(:, :, 4).*repmat(MESH.normal(1, :), nbound, 1) + sigma(:, :, 2).*repmat(MESH.normal(2, :), nbound, 1) + sigma(:, :, 6).*repmat(MESH.normal(3, :), nbound, 1);
TT(:,:,3)=sigma(:, :, 5).*repmat(MESH.normal(1, :), nbound, 1) + sigma(:, :, 6).*repmat(MESH.normal(2, :), nbound, 1) + sigma(:, :, 3).*repmat(MESH.normal(3, :), nbound, 1);

txty_el=zeros(nbound,ndim);
cont_integ=zeros(nbound,nn,ndim);
for i=1:nn
    txty_el(:,1)=TT(:,i,1);
    txty_el(:,2)=TT(:,i,2);
    txty_el(:,3)=TT(:,i,3);
    txty_int=N'*txty_el;
    
    NODES_el=NODES(:,SEGFA(:,i));
    NODES_temp=NODES_el-repmat(NODES_el(:,1),1,nbound);
    NODES_loc=zeros(2, nbound);
    l1=sqrt(NODES_temp(1,2)^2+NODES_temp(2,2)^2+NODES_temp(3,2)^2);
    l2=sqrt(NODES_temp(1,3)^2+NODES_temp(2,3)^2+NODES_temp(3,3)^2);
    theta=acos(dot(NODES_temp(:,2), NODES_temp(:,3))/(l1*l2));
    NODES_loc(1, 2)=l1;
    NODES_loc(2, 3)=l2*sin(theta);
    if nbound>3
        NODES_loc(1,4)=0.5*NODES_loc(1, 2);
        NODES_loc(2,4)=0.5*NODES_loc(2, 3);
        NODES_loc(2,5)=0.5*NODES_loc(2, 3);
        NODES_loc(1,6)=0.5*NODES_loc(1, 2);
        if nbound>6
            NODES_loc(1,7)=1/3*NODES_loc(1, 2);
            NODES_loc(2,7)=1/3*NODES_loc(2, 3);
        end
    end
    
    dNdU_ip=dNdU(:,:,1);
    J=NODES_loc*dNdU_ip;
    detJ=det(J);
    
    cont_integ(:,i,1)=N*(ipw(:).*txty_int(:,1))*detJ;
    cont_integ(:,i,2)=N*(ipw(:).*txty_int(:,2))*detJ;
    cont_integ(:,i,3)=N*(ipw(:).*txty_int(:,3))*detJ;
    
end

cont_integ_x=cont_integ(:,:,1);
cont_integ_y=cont_integ(:,:,2);
cont_integ_z=cont_integ(:,:,3);
cont_tx=accumarray(SEGFA(:), cont_integ_x(:), [nnod,1]);
cont_ty=accumarray(SEGFA(:), cont_integ_y(:), [nnod,1]);
cont_tz=accumarray(SEGFA(:), cont_integ_z(:), [nnod,1]);


F=zeros(ndim*nnod, 1);
F(1:ndim:(end-(ndim-1)))=cont_tx;
F(2:ndim:(end-(ndim-2)))=cont_ty;
F(3:ndim:end)=cont_tz;



