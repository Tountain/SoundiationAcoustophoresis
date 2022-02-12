function snm_ib = scattering_coeff_irregular_body(BC, db_bs_coeff)
%%
% this function calculate the scattering coefficient of irregular body
% "snm_ib" under Hard wall (Neumann) Boundary Conditions.
% For irregular body, the scattering coefficients are related to the
% beam-shape coefficient 'db_bs_coeff'.
% The irregular body maps to spherical body using mapping coefficients
% 'Cn'.
%%

% [Gam, Gam_partial_u, Amd, Amd_partial_u, Rau, Rau_partial_u] = structual_functions(Cn, n, m);
% 
% snm_ib = (Rau * Gam_partial_u - Rau_partial_u * Gam) ...
%         / (Rau_partial_u * Amd - Rau * Amd_partial_u);

%%

% The saving form of 'db_bs_coeff' in "database_beam_scattering_coeffs.m"
% are:
%    db_bs_coeff = [a(0, 0)      0       0       0       0      0   ...
%                   a(1,-1) a(1, 0) a(1, 1)      0       0      0   ...
%                   a(2,-2) a(2,-1) a(2, 0) a(2, 1) a(2, 2)     0
%                     ...     ...     ...     ...     ...      ...  ...],
% which is a 'side-trangular' matrix.
%
% In this code, in order to simplify the coding and make the matrix
% straightforward, we temperally re-arrange the 'db_bs_coeff' into a
% 'central-trangular' matrix as:
%    bs_coeff_ctl = [0      ...     0         0      a(0, 0)      0         0       ...      0 
%                    0      ...     0      a(1,-1)   a(1, 0)   a(1, 1)      0       ...      0
%                    0      ...  a(2,-2)   a(2,-1)   a(2, 0)   a(2, 1)   a(2, 2)    ...      0
%                    ...    ...     ...       ...       ...       ...       ...     ...     ...].
%
% In this way, the 'bs_coeff_ctl' can be general expressed as:
%    bs_coeff_ctl = [a(0,-N)      ...      a(0,-1)   a(0, 0)   a(0, 1)       ...      a(0, N)
%                    a(1,-N)      ...      a(1,-1)   a(1, 0)   a(1, 1)       ...      a(1, N)
%                    a(2,-N)      ...      a(2,-1)   a(2, 0)   a(2, 1)       ...      a(2, N)
%                      ...        ...        ...       ...       ...         ...        ...
%                    a(N,-N)      ...      a(N,-1)   a(N, 0)   a(N, 1)       ...      a(N, N)].
% NOTE that if the first index smaller than the absolution of the second
% index, the corresponding elements are zeros, such a(1,-2), a(1, 2), and
% a(N-1, N).

[row, col] = size(db_bs_coeff);
temp = [];
N = row-1;
% extract the meaningful information from 'db_bs_coeff' and temperally
% save in 'temp'
for nn = 0 : N
    for mm = -nn : nn
        temp = [temp, db_bs_coeff(nn+1, nn+mm+1)];
    end
end

ii = 1;
jj = 1;
bs_coeff_ctl = zeros(1, row*col);
% rearrange the 'temp' to 'central-trangular' matrix (if nn < abs(mm), the corresponding data is zero)
for nn = 0 : N
    for mm = -N : N
        if nn >= abs(mm)
            bs_coeff_ctl(jj) = temp(ii);
            ii = ii + 1;
            jj = jj + 1;
        else
            jj = jj + 1;
        end
    end
end
bs_coeff_ctl = reshape(bs_coeff_ctl, col, row).';


%% database of system prepare

A = zeros(row, col);        % eigen vector
Amd_database = cell(row, col);
% % ================== Single thread version ==================
% for n_s = 0 : N
%     for m_s = -N : N
%         [Gam, Gam_partial_u, Amd, Amd_partial_u, ~, ~] = structual_functions(n_s, m_s);
%         if strcmp(BC, 'rigid') == 1
%             A(n_s+1, m_s+N+1) = (bs_coeff_ctl(:, m_s+N+1)).' * Gam_partial_u;           % summation of index 'nn'
%             Amd_database{n_s+1, m_s+N+1} = Amd_partial_u;
%         else
%             A(n_s+1, m_s+N+1) = (bs_coeff_ctl(:, m_s+N+1)).' * Gam;                     % summation of index 'nn'
%             Amd_database{n_s+1, m_s+N+1} = Amd;
%         end
%     end
%     fprintf('Irregular object Scattering Coefficients Database Preparing %d%% \n',round(100*n_s/N));
% end
% ===================== Parallel version 1 ===================
parfor n_s = 0 : N
    %indices_1 = n_s+1;
    temp_A = zeros(1, col);
    temp_Amd_database = cell(1, col);
    for m_s = -N : N
        %indices_2 = N+1;
        [Gam, Gam_partial_u, Amd, Amd_partial_u, ~, ~] = structual_functions(n_s, m_s);
        if strcmp(BC, 'rigid') == 1
            temp_A(1, m_s+N+1) = (bs_coeff_ctl(:, m_s+N+1)).' * Gam_partial_u;           % summation of index 'nn'
            temp_Amd_database{1, m_s+N+1} = Amd_partial_u;
        else
            temp_A(1, m_s+N+1) = (bs_coeff_ctl(:, m_s+N+1)).' * Gam;                     % summation of index 'nn'
            temp_Amd_database{1, m_s+N+1} = Amd;
        end
    end
    A(n_s+1, :) = temp_A;
    Amd_database(n_s+1, :) = temp_Amd_database;
    fprintf('Irregular object Scattering Coefficients Database Preparing %d%% \n',round(100*n_s/N));
end
% % ===================== Parallel version 2 ===================
% for n_s = 0 : N
%     %indices_1 = n_s+1;
%     temp_A = zeros(1, col);
%     temp_Amd_database = cell(1, col);
%     indices_2 = N+1;
%     parfor m_s = -N : N
%         [Gam, Gam_partial_u, Amd, Amd_partial_u, ~, ~] = structual_functions(n_s, m_s);
%         if strcmp(BC, 'rigid') == 1
%             temp_A(1, m_s+indices_2) = (bs_coeff_ctl(:, m_s+indices_2)).' * Gam_partial_u;           % summation of index 'nn'
%             temp_Amd_database{1, m_s+indices_2} = Amd_partial_u;
%         else
%             temp_A(1, m_s+indices_2) = (bs_coeff_ctl(:, m_s+indices_2)).' * Gam;                     % summation of index 'nn'
%             temp_Amd_database{1, m_s+indices_2} = Amd;
%         end
%     end
%     A(n_s+1, :) = temp_A;
%     Amd_database(n_s+1, :) = temp_Amd_database;
%     fprintf('Irregular object Scattering Coefficients Database Preparing %d%% \n',round(100*n_s/N));
% end

% all eigen matrix
Mat = cell(2*N+1,1);
for m_s = -N : N
    for n_s = 0 : N
        %(bs_coeff_ctl(:, m_s+N+1)).' .* (Amd_database{n_s+1, m_s+N+1}).';
    	Mat{m_s+N+1} = [Mat{m_s+N+1}; (bs_coeff_ctl(:, m_s+N+1)).' .* (Amd_database{n_s+1, m_s+N+1}).'];
    end
end


%% solving the system equations

snm_ib_ctl = zeros(row, col);
for m_s = -N : N
    snm_ib_ctl(:, m_s+N+1) = - pinv(Mat{m_s+N+1}) * A(:, m_s+N+1);
%     snm_ib_ctl(:, m_s+N+1) = - (Mat{m_s+N+1}) \ A(:, m_s+N+1);
end

%% transfer the 'central-trangular' matrix to 'side-trangular' matrix of scattering coefficients
% this transfer is essential as the later calculations are all based on
% 'side-trangular' form, which can save the computational resource.
jj = 1;
temp = [];
snm_ib_ctl = reshape(snm_ib_ctl.', 1, row*col);
for nn = 0 : N
    for mm = -N : N
        if nn >= abs(mm)
            temp = [temp snm_ib_ctl(jj)];
            jj = jj + 1;
        else
            jj = jj + 1;
        end
    end
end 

ii = 1;
snm_ib = zeros(row, col);
for nn = 0 : N
    for mm = -nn : nn
        snm_ib(nn+1, nn+mm+1) = temp(ii);
        ii = ii + 1;
    end
end

%%
