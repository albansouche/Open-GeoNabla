function [A_all, InvJ_all, DetJ_all] = matrix_entries_EP(MESH, D_el_ip, PARAM, SOLVER, elast)

% Part of Open-GeoNabla, copyright GPLv3, 2018
% https://github.com/albansouche/Open-GeoNabla/
% Physics of Geological Processes (PGP) , The NJORD Centre, Dept of Geosciences, University of Oslo
% Author: Alban Souche


% 3D VERSION OF THE MILAMIN MECHANICAL2D.M:
% %   MECHANICAL2D Two dimensional finite element mechanical problem solver of MILAMIN
% %   Part of MILAMIN: MATLAB-based FEM solver for large problems
% %   Version 1.0.1
% %   Copyright (C) 2011, M. Dabrowski, M. Krotkiewski, D.W. Schmid
% %   University of Oslo, Physics of Geological Processes
% %   http://milamin.org
% %   See License file for terms of use.



%==========================================================================
% MODEL INFO
%==========================================================================
[ndim,     ~] = size(MESH.NODES);
[nnodel, nel] = size(MESH.ELEMS);
nedof         = nnodel*ndim;
nip           = SOLVER.nip;

%==========================================================================
% BLOCKING PARAMETERS (nelblo must be < nel)
%==========================================================================
nelblo          = min(nel, SOLVER.nelblo);
nblo            = ceil(nel/nelblo);

%==========================================================================
% i) PREPARE INTEGRATION POINTS & DERIVATIVES wrt LOCAL COORDINATES
%==========================================================================
[IP_X, IP_w]    = ip_tetra(nip);
[   ~, dNdu]    = shp_deriv_tetra(IP_X, nnodel);

%==================================================================
% DECLARE VARIABLES (ALLOCATE MEMORY)
%==================================================================
if or(elast==1, isequal(PARAM.ep_phi,PARAM.ep_xsi))
    A_block     = zeros(nelblo, nedof*(nedof+1)/2);
    A_all       = zeros(nedof*(nedof+1)/2,nel);
elseif max(PARAM.ep_phi)>min(PARAM.ep_xsi)
    A_block     = zeros(nelblo, nedof*nedof);
    A_all       = zeros(nedof*nedof,nel);
end
InvJ_all    = zeros(nel,ndim*ndim);
DetJ_all    = zeros(nel,1);


il          = 1;
iu          = nelblo;
%==================================================================
% i) BLOCK LOOP - MATRIX COMPUTATION
%==================================================================

for ib = 1:nblo
    %==============================================================
    % ii) FETCH DATA OF ELEMENTS IN BLOCK
    %==============================================================
    ECOORD_x = reshape( MESH.NODES(1,MESH.ELEMS(:,il:iu)), nnodel, nelblo);
    ECOORD_y = reshape( MESH.NODES(2,MESH.ELEMS(:,il:iu)), nnodel, nelblo);
    ECOORD_z = reshape( MESH.NODES(3,MESH.ELEMS(:,il:iu)), nnodel, nelblo);

    %==============================================================
    % iii) INTEGRATION LOOP
    %==============================================================
    A_block(:)    = 0;
    

    for ip=1:nip

        %==========================================================
        % iv) LOAD SHAPE FUNCTIONS DERIVATIVES FOR INTEGRATION POINT
        %==========================================================
%         Ni      =        N{ip};
        dNdui   =     dNdu{ip};
        
        %==========================================================
        % v) CALCULATE JACOBIAN, ITS DETERMINANT AND INVERSE
        %==========================================================
        Jx          = ECOORD_x'*dNdui;
        Jy          = ECOORD_y'*dNdui;
        Jz          = ECOORD_z'*dNdui;
       
        detJ        = Jx(:,1).*Jy(:,2).*Jz(:,3) + ...
                      Jx(:,2).*Jy(:,3).*Jz(:,1) + ...
                      Jx(:,3).*Jy(:,1).*Jz(:,2) - ...
                      Jx(:,3).*Jy(:,2).*Jz(:,1) - ...
                      Jx(:,2).*Jy(:,1).*Jz(:,3) - ...
                      Jx(:,1).*Jy(:,3).*Jz(:,2) ;
        
        invdetJ     = 1.0./detJ;
        det_aa = Jy(:,2).*Jz(:,3) - Jz(:,2).*Jy(:,3);
        invJx(:,1)  = det_aa.*invdetJ;
        det_ab = Jx(:,3).*Jz(:,2) - Jz(:,3).*Jx(:,2);
        invJy(:,1)  = det_ab.*invdetJ;
        det_ac = Jx(:,2).*Jy(:,3) - Jy(:,2).*Jx(:,3);
        invJz(:,1)  = +det_ac.*invdetJ;
        det_ba = Jy(:,3).*Jz(:,1) - Jz(:,3).*Jy(:,1);
        invJx(:,2)  = +det_ba.*invdetJ;
        det_bb = Jx(:,1).*Jz(:,3) - Jz(:,1).*Jx(:,3);
        invJy(:,2)  = +det_bb.*invdetJ;
        det_bc = Jx(:,3).*Jy(:,1) - Jy(:,3).*Jx(:,1);
        invJz(:,2)  = +det_bc.*invdetJ;
        det_ca = Jy(:,1).*Jz(:,2) - Jz(:,1).*Jy(:,2);
        invJx(:,3)  = +det_ca.*invdetJ;
        det_cb = Jx(:,2).*Jz(:,1) - Jz(:,2).*Jx(:,1);
        invJy(:,3)  = +det_cb.*invdetJ;
        det_cc = Jx(:,1).*Jy(:,2) - Jy(:,1).*Jx(:,2);
        invJz(:,3)  = +det_cc.*invdetJ;
       
        %==========================================================
        % vi) DERIVATIVES wrt GLOBAL COORDINATES
        %==========================================================
        dNdx        = invJx*dNdui';
        dNdy        = invJy*dNdui';
        dNdz        = invJz*dNdui';
      
        %==========================================================
        % vii) NUMERICAL INTEGRATION OF ELEMENT MATRICES
        %==========================================================
        weight      = IP_w(ip)*detJ;
                    
        % ------------------------A matrix-------------------------
        if elast==1 %elatic
            
            % Load elasto-plastic tangent matrix entries

            D1 = D_el_ip(il:iu,ip,1);
            D2 = D_el_ip(il:iu,ip,2);
            D3 = D_el_ip(il:iu,ip,22);
            
            indx  = 1;
            for i = 1:nnodel
                % x-velocity equation
                for j = i:nnodel
                    A_block(:,indx) = A_block(:,indx) + ( D1.*dNdx(:,i).*dNdx(:,j) + D3.*dNdy(:,i).*dNdy(:,j) + D3.*dNdz(:,i).*dNdz(:,j)).*weight;
                    indx = indx+1;
                    A_block(:,indx) = A_block(:,indx) + ( D2.*dNdx(:,i).*dNdy(:,j) + D3.*dNdy(:,i).*dNdx(:,j)).*weight;
                    indx = indx+1;
                    A_block(:,indx) = A_block(:,indx) + ( D2.*dNdx(:,i).*dNdz(:,j) + D3.*dNdz(:,i).*dNdx(:,j)).*weight;
                    indx = indx+1;
                end                
                % y-velocity equation
                for j = i:nnodel
                    if(j>i)
                        A_block(:,indx) = A_block(:,indx) + ( D2.*dNdy(:,i).*dNdx(:,j) + D3.*dNdx(:,i).*dNdy(:,j)).*weight;
                        indx = indx+1;
                    end
                    A_block(:,indx) = A_block(:,indx) + ( D1.*dNdy(:,i).*dNdy(:,j) + D3.*dNdx(:,i).*dNdx(:,j) + D3.*dNdz(:,i).*dNdz(:,j)).*weight;
                    indx = indx+1;
                    A_block(:,indx) = A_block(:,indx) + ( D2.*dNdy(:,i).*dNdz(:,j) + D3.*dNdz(:,i).*dNdy(:,j)).*weight;
                    indx = indx+1;
                end
                % z-velocity equation
                for j = i:nnodel
                    if(j>i)
                        A_block(:,indx) = A_block(:,indx) + ( D2.*dNdz(:,i).*dNdx(:,j) + D3.*dNdx(:,i).*dNdz(:,j)).*weight;
                        indx = indx+1;
                        A_block(:,indx) = A_block(:,indx) + ( D2.*dNdz(:,i).*dNdy(:,j) + D3.*dNdy(:,i).*dNdz(:,j)).*weight;
                        indx = indx+1;
                    end
                    A_block(:,indx) = A_block(:,indx) + ( D1.*dNdz(:,i).*dNdz(:,j) + D3.*dNdy(:,i).*dNdy(:,j) + D3.*dNdx(:,i).*dNdx(:,j)).*weight;
                    indx = indx+1;
                end
                
                
            end
            
        elseif PARAM.ep_phi==PARAM.ep_xsi % elastoplastic & associative plastic flow rule
            
            d11 = D_el_ip(il:iu,ip,1);
            d12 = D_el_ip(il:iu,ip,2);
            d13 = D_el_ip(il:iu,ip,3);
            d14 = D_el_ip(il:iu,ip,4);
            d15 = D_el_ip(il:iu,ip,5);
            d16 = D_el_ip(il:iu,ip,6);
            d22 = D_el_ip(il:iu,ip,8);
            d23 = D_el_ip(il:iu,ip,9);
            d24 = D_el_ip(il:iu,ip,10);
            d25 = D_el_ip(il:iu,ip,11);
            d26 = D_el_ip(il:iu,ip,12);
            d33 = D_el_ip(il:iu,ip,15);
            d34 = D_el_ip(il:iu,ip,16);
            d35 = D_el_ip(il:iu,ip,17);
            d36 = D_el_ip(il:iu,ip,18);
            d44 = D_el_ip(il:iu,ip,22);
            d45 = D_el_ip(il:iu,ip,23);
            d46 = D_el_ip(il:iu,ip,24);
            d55 = D_el_ip(il:iu,ip,29);
            d56 = D_el_ip(il:iu,ip,30);
            d66 = D_el_ip(il:iu,ip,36);
            
            indx  = 1;
            
            for i = 1:nnodel
               % x-velocity equation
                for j = i:nnodel
                    A_block(:,indx) = A_block(:,indx) + ( d11.*dNdx(:,i).*dNdx(:,j) + d14.*dNdy(:,i).*dNdx(:,j) + d16.*dNdz(:,i).*dNdx(:,j) ...
                                                        + d14.*dNdx(:,i).*dNdy(:,j) + d44.*dNdy(:,i).*dNdy(:,j) + d46.*dNdz(:,i).*dNdy(:,j) ...
                                                        + d16.*dNdx(:,i).*dNdz(:,j) + d46.*dNdy(:,i).*dNdz(:,j) + d66.*dNdz(:,i).*dNdz(:,j) ).*weight;
                    indx = indx+1;
                    A_block(:,indx) = A_block(:,indx) + ( d12.*dNdx(:,i).*dNdy(:,j) + d24.*dNdy(:,i).*dNdy(:,j) + d26.*dNdz(:,i).*dNdy(:,j) ...
                                                        + d14.*dNdx(:,i).*dNdx(:,j) + d44.*dNdy(:,i).*dNdx(:,j) + d46.*dNdz(:,i).*dNdx(:,j) ...
                                                        + d15.*dNdx(:,i).*dNdz(:,j) + d45.*dNdy(:,i).*dNdz(:,j) + d56.*dNdz(:,i).*dNdz(:,j) ).*weight;
                    indx = indx+1;
                    A_block(:,indx) = A_block(:,indx) + ( d13.*dNdx(:,i).*dNdz(:,j) + d34.*dNdy(:,i).*dNdz(:,j) + d36.*dNdz(:,i).*dNdz(:,j) ...
                                                        + d15.*dNdx(:,i).*dNdy(:,j) + d45.*dNdy(:,i).*dNdy(:,j) + d56.*dNdz(:,i).*dNdy(:,j) ...
                                                        + d16.*dNdx(:,i).*dNdx(:,j) + d46.*dNdy(:,i).*dNdx(:,j) + d66.*dNdz(:,i).*dNdx(:,j) ).*weight;
                    indx = indx+1;
                end                
                % y-velocity equation
                for j = i:nnodel
                    if(j>i)
                        A_block(:,indx) = A_block(:,indx) + ( d12.*dNdy(:,i).*dNdx(:,j) + d14.*dNdx(:,i).*dNdx(:,j) + d15.*dNdz(:,i).*dNdx(:,j) ...
                                                            + d24.*dNdy(:,i).*dNdy(:,j) + d44.*dNdx(:,i).*dNdy(:,j) + d45.*dNdz(:,i).*dNdy(:,j) ...
                                                            + d26.*dNdy(:,i).*dNdz(:,j) + d46.*dNdx(:,i).*dNdz(:,j) + d56.*dNdz(:,i).*dNdz(:,j) ).*weight;
                        indx = indx+1;
                    end
                    A_block(:,indx) = A_block(:,indx) + ( d22.*dNdy(:,i).*dNdy(:,j) + d24.*dNdx(:,i).*dNdy(:,j) + d25.*dNdz(:,i).*dNdy(:,j) ...
                                                        + d24.*dNdy(:,i).*dNdx(:,j) + d44.*dNdx(:,i).*dNdx(:,j) + d45.*dNdz(:,i).*dNdx(:,j) ...
                                                        + d25.*dNdy(:,i).*dNdz(:,j) + d45.*dNdx(:,i).*dNdz(:,j) + d55.*dNdz(:,i).*dNdz(:,j) ).*weight;
                    indx = indx+1;
                    A_block(:,indx) = A_block(:,indx) + ( d23.*dNdy(:,i).*dNdz(:,j) + d34.*dNdx(:,i).*dNdz(:,j) + d35.*dNdz(:,i).*dNdz(:,j) ...
                                                        + d25.*dNdy(:,i).*dNdy(:,j) + d45.*dNdx(:,i).*dNdy(:,j) + d55.*dNdz(:,i).*dNdy(:,j) ...
                                                        + d26.*dNdy(:,i).*dNdx(:,j) + d46.*dNdx(:,i).*dNdx(:,j) + d56.*dNdz(:,i).*dNdx(:,j) ).*weight;
                    indx = indx+1;
                end
                % z-velocity equation
                for j = i:nnodel
                    if(j>i)
                        A_block(:,indx) = A_block(:,indx) + ( d13.*dNdz(:,i).*dNdx(:,j) + d15.*dNdy(:,i).*dNdx(:,j) + d16.*dNdx(:,i).*dNdx(:,j) ...
                                                            + d34.*dNdz(:,i).*dNdy(:,j) + d45.*dNdy(:,i).*dNdy(:,j) + d46.*dNdx(:,i).*dNdy(:,j) ...
                                                            + d36.*dNdz(:,i).*dNdz(:,j) + d56.*dNdy(:,i).*dNdz(:,j) + d66.*dNdx(:,i).*dNdz(:,j) ).*weight;
                        indx = indx+1;
                        A_block(:,indx) = A_block(:,indx) + ( d23.*dNdz(:,i).*dNdy(:,j) + d25.*dNdy(:,i).*dNdy(:,j) + d26.*dNdx(:,i).*dNdy(:,j) ...
                                                            + d34.*dNdz(:,i).*dNdx(:,j) + d45.*dNdy(:,i).*dNdx(:,j) + d46.*dNdx(:,i).*dNdx(:,j) ...
                                                            + d35.*dNdz(:,i).*dNdz(:,j) + d55.*dNdy(:,i).*dNdz(:,j) + d56.*dNdx(:,i).*dNdz(:,j) ).*weight;
                        indx = indx+1;
                    end
                    A_block(:,indx) = A_block(:,indx) + ( d33.*dNdz(:,i).*dNdz(:,j) + d35.*dNdy(:,i).*dNdz(:,j) + d36.*dNdx(:,i).*dNdz(:,j) ...
                                                        + d35.*dNdz(:,i).*dNdy(:,j) + d55.*dNdy(:,i).*dNdy(:,j) + d56.*dNdx(:,i).*dNdy(:,j) ...
                                                        + d36.*dNdz(:,i).*dNdx(:,j) + d56.*dNdy(:,i).*dNdx(:,j) + d66.*dNdx(:,i).*dNdx(:,j) ).*weight;
                    indx = indx+1;
                end
                
            end
            
        elseif PARAM.ep_phi>PARAM.ep_xsi 
            
            fprintf(1, 'Non associative plasticity not implemented. Use standard matrice assembly');

        end
        
    end
    
    %==============================================================
    % ix) WRITE DATA INTO GLOBAL STORAGE
    %==============================================================
    A_all(:,il:iu)    = A_block';
    InvJ_all(il:iu,:) = [invJx invJy invJz];
    DetJ_all(il:iu,:) = detJ;
    
    %==============================================================
    % READJUST START, END AND SIZE OF BLOCK. REALLOCATE MEMORY
    %==============================================================
    il  = il+nelblo;
    if(ib==nblo-1)
        nelblo      = nel-iu;
        if or(elast==1, isequal(PARAM.ep_phi,PARAM.ep_xsi))
            A_block     = zeros(nelblo, nedof*(nedof+1)/2);
        elseif max(PARAM.ep_phi)>min(PARAM.ep_xsi)
            A_block     = zeros(nelblo, nedof*nedof);
        end
        invJx       = zeros(nelblo, ndim);
        invJy       = zeros(nelblo, ndim);
        invJz       = zeros(nelblo, ndim);
        detJ        = zeros(nelblo, 1);
    end
    iu  = iu+nelblo;
end

   
