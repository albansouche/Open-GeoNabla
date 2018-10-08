function [F,cont_tx,cont_ty] = ContourTraction(MESH, SIGMA, seg_id, seg_id_corner)

% Part of Open-GeoNabla, copyright GPLv3, 2018
% https://github.com/albansouche/Open-GeoNabla/
% Physics of Geological Processes (PGP) , The NJORD Centre, Dept of Geosciences, University of Oslo
% Authors: Jan Cornet


% ContourTraction computes the contour integral of tractions around a rectangle
%    cont_integ=ContourTraction(MESH,SIGMA,nip)

%   *******************
%   *                 *    * nodes not considered for traction calculation
%   *                 *    # corner nodes shared
%   #                 *    . nodes considered for traction calculation
%   .                 *
%   .                 *
%   .                 *
%   .                 *
%   #                 *
%   *                 *
%   *                 *
%   *******************

% Traction only along specific segment number %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Correction from shared nodes at the corner
if and(seg_id,seg_id_corner)
    dummy1 = MESH.SEGMENTS(:,MESH.segment_markers==seg_id);
    dummy2 = MESH.SEGMENTS(:,MESH.segment_markers==seg_id_corner);
    nodes_shared_corner = dummy1(ismember(dummy1,dummy2));
    indx_seg_out = ~(MESH.segment_markers==seg_id);
    MESH.segment_markers(indx_seg_out) = [];
    MESH.SEGMENTS(:,indx_seg_out)      = [];
    MESH.normal(:,indx_seg_out)        = [];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NODES=MESH.NODES;
[ndim,nnod]=size(NODES);
nip=6;

SEGFA=MESH.SEGMENTS;
[nbound, nn]=size(SEGFA);
[ipuv, ipw] = ip_gauss(nip, 1);
[N,~] = shp_deriv_bar(ipuv,nbound);

sigma=zeros(nbound, nn, ndim*(ndim+1)/2);
sigma(:, :, 1)=reshape(SIGMA(1, SEGFA), size(SEGFA)); 
sigma(:, :, 2)=reshape(SIGMA(2, SEGFA), size(SEGFA)); 
sigma(:, :, 3)=reshape(SIGMA(3, SEGFA), size(SEGFA)); 

TT=zeros(nbound,nn,ndim);
%compute the components of the tractions
TT(:,:,1)=sigma(:, :, 1).*repmat(MESH.normal(1, :), nbound, 1) + sigma(:, :, 3).*repmat(MESH.normal(2, :), nbound, 1);
TT(:,:,2)=sigma(:, :, 3).*repmat(MESH.normal(1, :), nbound, 1) + sigma(:, :, 2).*repmat(MESH.normal(2, :), nbound, 1);

%perform the contour integral
txty_el=zeros(nbound,ndim);
cont_integ=zeros(nbound,nn,ndim);

for i=1:nn
    txty_el(:,1)=TT(:,i,1);
    txty_el(:,2)=TT(:,i,2);
    txty_int=N'*txty_el;
    
    detJ=abs(NODES(1,SEGFA(2,i))-NODES(1,SEGFA(1,i))...
        +1i*(NODES(2,SEGFA(2,i))-NODES(2,SEGFA(1,i))))/2;
    
    cont_integ(:,i,1)=N*(ipw(:).*txty_int(:,1))*detJ;
    cont_integ(:,i,2)=N*(ipw(:).*txty_int(:,2))*detJ;
end
   
cont_integ_x=cont_integ(:,:,1);
cont_integ_y=cont_integ(:,:,2);
cont_tx=accumarray(SEGFA(:), cont_integ_x(:), [nnod,1]);
cont_ty=accumarray(SEGFA(:), cont_integ_y(:), [nnod,1]);

F=zeros(ndim*nnod, 1);
F(1:ndim:end)=cont_tx;
F(2:ndim:end)=cont_ty;



% %%%%%%% FIX CORNER NODES FOR 2D %%%%%%%%
% if ndim==2
%     dof_corner_nodes = [ ndim*(nodes_shared_corner(:)-1)+1 ; ndim*(nodes_shared_corner(:)-1)+2 ];
%     F(dof_corner_nodes) = 2*F(dof_corner_nodes);
% end





    function [ipx, ipw] = ip_gauss1d(N)
        %  http://www.comlab.ox.ac.uk/people/nick.trefethen/gauss.m
        beta = .5./sqrt(1-(2*(1:N-1)).^(-2));
        T = diag(beta,1) + diag(beta,-1);
        [V,D] = eig(T);
        ipx = diag(D); [ipx,i] = sort(ipx);
        ipw = 2*V(1,i).^2;
    end

    function [ipx, ipw, siz] = ip_gauss(no_pts, ndim)
        %------------------------------------------------------
        %Gauss-Legendre quadrature - integration points and weights
        %------------------------------------------------------
        %initialize
        siz = no_pts^ndim;
        ipx=zeros(siz,ndim);
        ipw=zeros(siz,1);
        iter = 1;
        %use 1D GL
        [ipx1D,ipw1D]=ip_gauss1d(no_pts);
        ipx = ipx1D;
        ipw = ipw1D(:);
        
    end

    function [N, dNdu] = shp_deriv_bar(ipx, nnodel)
        nip  = size(ipx,1);
        N    = cell(nip,1);
        dNdu = cell(nip,1);
        for ii=1:nip
            x = ipx(ii);
            switch nnodel
                case 2
                    SHP = [(1-x)/2; (1+x)/2];
                    DERIV = [-1/2 1/2]; %wrt x
                case 3
                    SHP = [(x-1)/2*x; (1+x)/2*x; (1-x)*(x+1)];
                    DERIV = [x-1/2 x+1/2 -2*x]; %wrt x 
            end
            N{ii} = SHP;
            dNdu{ii} = DERIV';
        end
        %convert cells to arrays
        nip = size(ipx,1);
        N = cell2mat(N);
        N = reshape(N,nnodel,nip);
        dNdu = reshape(cell2mat(dNdu),nnodel,nip,1);
        dNdu = permute(dNdu, [1 3 2]);
    end

end