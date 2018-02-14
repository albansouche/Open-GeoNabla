function Disp = mechanical_solver(MESH, A_all, BC, rhs, SOLVER)

% Part of Open-GeoNabla, copyright GPLv3, 2018
% https://github.com/albansouche/Open-GeoNabla/
% Physics of Geological Processes (PGP) , The NJORD Centre, Dept of Geosciences, University of Oslo
% Author: Alban Souche

% Function build upon MILAMIN Version 1.0.1, Mutils 0.4.2


[ndim, nnod] = size(MESH.NODES);
sdof        = ndim*nnod;

switch SOLVER.solver_type
    
    case 'direct_nonsym'
        
    case 'direct_sym'
        
    case 'direct_sym_chol'
        
    case 'PCG'

        % Space assembly
        opts_sparse.symmetric  = 1;
        opts_sparse.n_node_dof =ndim;
        opts_sparse.nthreads   = SOLVER.nb_cpu;
        opts_sparse.cpu_affinity = 0;
        A = sparse_create(MESH.ELEMS, A_all, opts_sparse);
        clear A_all;
        
        % BC
        Disp         = zeros(sdof  , 1);
        Disp(BC.ind) = BC.val;
        Free        = 1:sdof;
        Free(BC.ind)= [];
        Rhs         = rhs - A*Disp - A'*Disp;
        Rhs         = Rhs(Free);
        A           = A(Free, Free);
    
        % Preconditioning
        precond = 1./spdiags(A,0);
        A_conv = sparse_convert(A, opts_sparse);
        Ax = @(x) spmv(A_conv,x);

%         tic
        [Disp(Free),~,~, iter] = pcg(Ax, Rhs, SOLVER.ep_rel_res_tol, 10000, @(x) precond.*x);
%         fprintf(1, ['PCG(',num2str(iter),'it)    SOLVED in ',num2str(toc), ' sec. \n']);

   case 'AGMG'

   
end





