function [Edot, W]= strain_rate_calc(MESH, INVJ, V, pts_eva_loc, nelblo)

% Part of Open-GeoNabla, copyright GPLv3, 2018
% https://github.com/albansouche/Open-GeoNabla/
% Physics of Geological Processes (PGP) , The NJORD Centre, Dept of Geosciences, University of Oslo
% Author: Alban Souche


% Parameters
ndim         = size(MESH.NODES,1);
[nnodel,nel] = size(MESH.ELEMS);
npts         = length(pts_eva_loc);
nelblo       = min(nel, nelblo);
nblo         = ceil(nel/nelblo);
il           = 1;
iu           = nelblo;

% Derivative of shape functions
[ ~, dNdu]   = shp_deriv_triangle(pts_eva_loc, nnodel);

% Allocation
Edot.xx = zeros(npts,nel);
Edot.yy = zeros(npts,nel);
Edot.xy = zeros(npts,nel);
W.xy    = zeros(npts,nel);

for ib = 1:nblo
    
    % Extract inverse of Jacobian
    invJx = INVJ(il:iu,1:2);
    invJy = INVJ(il:iu,3:4);
    
    % Extract velocities
    Vx = reshape( V(ndim*(MESH.ELEMS(:,il:iu)-1)+1), nnodel,nelblo)';
    Vy = reshape( V(ndim*(MESH.ELEMS(:,il:iu)-1)+2), nnodel,nelblo)';
    
    for ieva=1:npts
        
        %Extract spatial derivative
        dNdui = dNdu{ieva};
        dNdx  = invJx*dNdui';
        dNdy  = invJy*dNdui';
        
        %Calculate strain rates (remember that deviatoric quantity is: Edot_dev(ij) = Edot(ij) - (1/3)*kro(ij)*Edot(kk) )
        Edot.xx(ieva,il:iu)  = 2/3*sum(dNdx.*Vx,2) - 1/3*sum(dNdy.*Vy,2);
        Edot.yy(ieva,il:iu)  = 2/3*sum(dNdy.*Vy,2) - 1/3*sum(dNdx.*Vx,2);
        Edot.xy(ieva,il:iu)  = 0.5*(sum(dNdy.*Vx,2)+sum(dNdx.*Vy,2));
        
        %Calculate rotational
        W.xy(ieva,il:iu)     = 0.5*(sum(dNdy.*Vx,2)-sum(dNdx.*Vy,2));
        
    end
    il  = il+nelblo;
    if(ib==nblo-1)
        nelblo 	  = nel-iu;
    end
    iu  = iu+nelblo;
end

% % add the numerical error from incompressiblity to the Edotzz (=0)
% Edotzz = (Edot.xx - Edot.yy);
Edot.zz = - (Edot.xx + Edot.yy);
Edot.nd2 = sqrt( (Edot.xx.^2 + Edot.yy.^2 + Edot.zz.^2)/2 + Edot.xy.^2);
% Edot.nd2 = sqrt( (Edot.xx.^2 + Edot.yy.^2)/2 + Edot.xy.^2);




