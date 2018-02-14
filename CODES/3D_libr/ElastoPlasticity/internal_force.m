function Fint = internal_force(MESH, INVJ, DetJ, Sdev, P, SOLVER) 

% Part of Open-GeoNabla, copyright GPLv3, 2018
% https://github.com/albansouche/Open-GeoNabla/
% Physics of Geological Processes (PGP) , The NJORD Centre, Dept of Geosciences, University of Oslo
% Author: Alban Souche


% Parameters
ndim         = size(MESH.NODES,1);
[nnodel,nel] = size(MESH.ELEMS);
nedof        = nnodel*ndim;
nip          = SOLVER.nip;
nelblo       = min(nel, SOLVER.nelblo);
nblo         = ceil(nel/nelblo);
il           = 1;
iu           = nelblo;

% Derivative of shape functions
[ipuv, ipw]  = ip_tetra(nip);
[ ~, dNdu]   = shp_deriv_tetra(ipuv, nnodel);

% Allocation
F_all   = zeros(nedof,nel); 
F_block = zeros(nelblo, nedof);

for ib = 1:nblo
    
    % Extract determinant of Jacobian
    detJ  = DetJ(il:iu);
    % Extract inverse of Jacobian
    invJx = INVJ(il:iu,1:3);
    invJy = INVJ(il:iu,4:6);
    invJz = INVJ(il:iu,7:9);
    % resert block vectors to zeros
    F_block(:) = 0;

    for ip=1:nip
        
        % Stresses at itegration points
        S1 = Sdev.xx(ip,il:iu)' + P(ip,il:iu)';
        S2 = Sdev.yy(ip,il:iu)' + P(ip,il:iu)';
        S3 = Sdev.zz(ip,il:iu)' + P(ip,il:iu)';
        S4 = Sdev.xy(ip,il:iu)';
        S5 = Sdev.yz(ip,il:iu)';
        S6 = Sdev.zx(ip,il:iu)';
        
        % Extract spatial derivative
        dNdui  = dNdu{ip};
        dNdx   = invJx*dNdui';
        dNdy   = invJy*dNdui';
        dNdz   = invJz*dNdui';
        
        % Weight of the itegration point
        weight = ipw(ip)*detJ;
       
        % Internal force vector
        F_block(:,1:ndim:nedof) =  F_block(:,1:ndim:nedof) + bsxfun(@times,(S1.*weight),dNdx) + bsxfun(@times,(S4.*weight),dNdy) + bsxfun(@times,(S6.*weight),dNdz);
        F_block(:,2:ndim:nedof) =  F_block(:,2:ndim:nedof) + bsxfun(@times,(S2.*weight),dNdy) + bsxfun(@times,(S4.*weight),dNdx) + bsxfun(@times,(S5.*weight),dNdz);
        F_block(:,3:ndim:nedof) =  F_block(:,3:ndim:nedof) + bsxfun(@times,(S3.*weight),dNdz) + bsxfun(@times,(S5.*weight),dNdy) + bsxfun(@times,(S6.*weight),dNdx);
        
    end
     
    F_all(:,il:iu)    = F_block';
    
    il  = il+nelblo;
    
    if(ib==nblo-1)
        nelblo 	= nel-iu;
        F_block = zeros(nelblo, nedof);
    end
    iu  = iu+nelblo;
end

ELEM_DOF = zeros(nedof, nel,'int32');
ELEM_DOF(1:ndim:end,:) = ndim*(MESH.ELEMS-1)+1;
ELEM_DOF(2:ndim:end,:) = ndim*(MESH.ELEMS-1)+2;
ELEM_DOF(3:ndim:end,:) = ndim*(MESH.ELEMS-1)+3;
Fint  = accumarray(ELEM_DOF(:), F_all(:));

