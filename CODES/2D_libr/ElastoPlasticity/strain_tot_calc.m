function Et = strain_tot_calc(MESH, INVJ, U, pts_eva_loc, nelblo)

% Part of Open-GeoNabla, copyright GPLv3, 2018
% https://github.com/albansouche/Open-GeoNabla/
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
[ ~, dNdu]   = shp_deriv_triangle(pts_eva_loc, nnodel);

% Allocation
Et.xx = zeros(npts,nel);
Et.yy = zeros(npts,nel);
Et.zz = zeros(npts,nel);
Et.xy = zeros(npts,nel);

for ib = 1:nblo
    
    % Extract inverse of Jacobian
    invJx = INVJ(il:iu,1:2);
    invJy = INVJ(il:iu,3:4);
    
    % Extract velocities
    Ux = reshape( U(ndim*(MESH.ELEMS(:,il:iu)-1)+1), nnodel,nelblo)';
    Uy = reshape( U(ndim*(MESH.ELEMS(:,il:iu)-1)+2), nnodel,nelblo)';
    
    for ieva=1:npts
        
        %Extract spatial derivative
        dNdui = dNdu{ieva};
        dNdx  = invJx*dNdui';
        dNdy  = invJy*dNdui';
        
        %Calculate strain rates (remember that deviatoric quantity is: Edot_dev(ij) = Edot(ij) - (1/3)*kro(ij)*Edot(kk) )
        Et.xx(ieva,il:iu) = sum(dNdx.*Ux,2);
        Et.yy(ieva,il:iu) = sum(dNdy.*Uy,2);
        Et.xy(ieva,il:iu) = 0.5*(sum(dNdy.*Ux,2)+sum(dNdx.*Uy,2));
        
    end
    il  = il+nelblo;
    if(ib==nblo-1)
        nelblo 	  = nel-iu;
    end
    iu  = iu+nelblo;
end




