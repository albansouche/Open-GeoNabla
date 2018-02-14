function Et = strain_tot_calc(MESH, INVJ, U, pts_eva_loc, nelblo)

% Part of OpenGeoNabla, copyright GPLv3, 2018
% https://github.com/albansouche/OpenGeoNabla/
% Physics of Geological Processes (PGP) , The NJORD Centre, Dept of Geosciences, University of Oslo
% Author: Alban Souche

% Parameters
ndim         = size(MESH.NODES,1);
[nnodel,nel] = size(MESH.ELEMS);
npts         = size(pts_eva_loc,1);
nelblo       = min(nel, nelblo);
nblo         = ceil(nel/nelblo);
il           = 1;
iu           = nelblo;

% Derivative of shape functions
[ ~, dNdu]   = shp_deriv_tetra(pts_eva_loc, nnodel);
 
% Allocation
Et.xx = zeros(npts,nel);
Et.yy = zeros(npts,nel);
Et.zz = zeros(npts,nel);
Et.xy = zeros(npts,nel);
Et.yz = zeros(npts,nel);
Et.zx = zeros(npts,nel);

for ib = 1:nblo
    
    % Extract inverse of Jacobian
    invJx = INVJ(il:iu,1:3);
    invJy = INVJ(il:iu,4:6);
    invJz = INVJ(il:iu,7:9);
    
    % Extract velocities
    Ux = reshape( U(ndim*(MESH.ELEMS(:,il:iu)-1)+1), nnodel,nelblo)';
    Uy = reshape( U(ndim*(MESH.ELEMS(:,il:iu)-1)+2), nnodel,nelblo)';
    Uz = reshape( U(ndim*(MESH.ELEMS(:,il:iu)-1)+3), nnodel,nelblo)';
    
    for ieva=1:npts
        
        % Extract derivative
        dNdui = dNdu{ieva};
        dNdx  = invJx*dNdui';
        dNdy  = invJy*dNdui';
        dNdz  = invJz*dNdui';
        
        % Calculate strain rates
        Et.xx(ieva,il:iu) = sum(dNdx.*Ux,2);
        Et.yy(ieva,il:iu) = sum(dNdy.*Uy,2);
        Et.zz(ieva,il:iu) = sum(dNdz.*Uz,2);
        Et.xy(ieva,il:iu) = 0.5*(sum(dNdy.*Ux,2)+sum(dNdx.*Uy,2));
        Et.yz(ieva,il:iu) = 0.5*(sum(dNdz.*Uy,2)+sum(dNdy.*Uz,2));
        Et.zx(ieva,il:iu) = 0.5*(sum(dNdx.*Uz,2)+sum(dNdz.*Ux,2));
        
    end
    il  = il+nelblo;
    if(ib==nblo-1)
        nelblo 	  = nel-iu;
    end
    iu  = iu+nelblo;
end




