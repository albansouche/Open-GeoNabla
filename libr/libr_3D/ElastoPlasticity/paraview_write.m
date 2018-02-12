function paraview_write(etype, NODES, ELEMS, PHASE, DATA, filename, itime)


%multiple DATAution arrays should be supported

if nargin<6
    itime = 0;
end

%--------------------------------------------------------------------------
%GENERAL
%--------------------------------------------------------------------------
otype = 'float';

%--------------------------------------------------------------------------
%MODEL INFO
%--------------------------------------------------------------------------
[ndim,  nnod] = size(NODES);
[nnodel, nel] = size(ELEMS);

% % etype = element_type(ndim, nnodel);
% if ndim==2
%    if nnodel==3 ; etype = 'tri3' ; end
%    if nnodel==6 ; etype = 'tri6' ; end
% elseif ndim==3
%    if nnodel==4 ; etype = 'tet4' ; end 
%    if nnodel==10; etype = 'tet10'; end 
%    if nnodel==8 ; etype = 'quad4'; end    
% end

%--------------------------------------------------------------------------
%INPUT VALIDATION
%--------------------------------------------------------------------------
if any(strcmp(etype,{'tri7', 'quad9','tet15','brick27'}))
    error('Element not supported.')
    %we could introduce here
    %1) element conversion (tet15->tet10, brick27->brick20, quad9->quad8)
    %2) element splitting (e.g., quad9->quad4)
    %3) data size reduction (only corner nodes)
end

nd = size(DATA,2);
switch nd
    case nnod
        data_type = 'point';
    case nel
        data_type = 'cell';
    otherwise
        error('Incorrect data size.')
end

%--------------------------------------------------------------------------
%ADD THIRD DIMENSION IF NEEDED
%--------------------------------------------------------------------------
if ndim == 2
    NODES = [NODES; zeros(1,nnod)];
    if size(DATA,1) == 2
        DATA = [DATA; zeros(1,size(DATA,2))];
    end
end

%--------------------------------------------------------------------------
%ADJUST CONVENTIONS
%--------------------------------------------------------------------------
%FIX OTHER CASES
switch etype
    case 'tri6'
    ELEMS(4:6,:) = ELEMS([6 4 5],:);
end

filename = strcat(filename,'_');
fid = fopen(strcat(strcat(filename,num2str(itime, '%.4d')),'.vtk'),'w','b');
fprintf(fid,'# vtk DataFile Version 3.0\n');
fprintf(fid,'my cool data\n');
fprintf(fid,'binary\n');
fprintf(fid,'DATASET UNSTRUCTURED_GRID\n');

%NODES
fprintf(fid,['POINTS %d ',otype,'\n'], nnod);
fwrite(fid, NODES, otype);

%ELEMS
%we should add support for points and lines
fprintf(fid,'\nCELLS %d %d\n', nel, nel*(nnodel+1));
ELEMS = ELEMS-1;
ELEMS = [nnodel*ones(1,nel); ELEMS];
fwrite(fid, ELEMS, 'int');

fprintf(fid,'\nCELL_TYPES %d\n', nel);
switch etype
    case 'tri3'
        cell_type =  5;
    case 'tri6'
        cell_type = 22;
    case 'quad4'
        cell_type =  9;
    case 'quad8'
        cell_type = 23;
    case 'tet4'
        cell_type = 10;
    case 'tet10'
        cell_type = 24;
    case 'brick8'
        cell_type = 12;
    case 'brick20'        
        cell_type = 25;
end

fwrite(fid, cell_type*ones(nel,1), 'int');

%PHASE
if ~isempty(PHASE)
    fprintf(fid,'\nCELL_DATA %d\n', nel);
    fprintf(fid,'SCALARS phase int 1\n');
    fprintf(fid,'LOOKUP_TABLE default\n');
    fwrite(fid, PHASE, 'int');
end

%DATA
switch data_type
    case 'point'        
        fprintf(fid,'\nPOINT_DATA %d\n', nnod);        
    case 'cell'        
        fprintf(fid,'\nCELL_DATA %d\n', nel); %fix jesli jest juz cell_data (phases) to nie moze sie powtarzac etykieta cell_data       
end
if size(DATA,1)==1
    fprintf(fid,['SCALARS result ',otype,' 1\n']);
    fprintf(fid,'LOOKUP_TABLE default\n');
elseif size(DATA,1)==3
    fprintf(fid,['VECTORS result ',otype,'\n']);
else
    fprintf(fid,['TENSORS result ',otype,'\n']);    
end
fwrite(fid, DATA, otype);
fclose(fid);
