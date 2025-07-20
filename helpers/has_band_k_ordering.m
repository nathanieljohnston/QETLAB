%%  HAS_BAND_K_ORDERING    Determines whether a matrix has bandwidth ≤ k up to symmetric permutation
%   This function has two required arguments:
%     A: a structurally symmetric matrix
%     k: a non-negative integer
%
%   This function has two output arguments:
%     has_ordering: whether there exists an ordering of A with bandwidth ≤ k
%     ordering: the ordering, if one exists; otherwise, just an empty vector
%
%   [has_ordering, ordering] = has_band_k_ordering(A, k) returns has_ordering,
%   ordering such that if has_ordering is true, then A(ordering, ordering) has
%   bandwidth ≤ k; otherwise, ordering is an empty vector. Recall that the
%   bandwidth of an n×n matrix A is the minimum non-negative integer k ∈ [0, n]
%   such that A(i, j) = 0 whenever |i - j| ≥ k (thus, zero matrices have
%   bandwidth 0, diagonal matrices have bandwidth 1, and so on). This function,
%   therefore, determines whether A has bandwidth at most k up to symmetric
%   permutation (reordering of rows and columns in the same way).
%
%   This function has no optional arguments.
%
%   The original minimization algorithm (which calls the recognition procedure
%   defined here for incrementing values of k until the true minimum is found)
%   is presented in Section 2 of [1]. The lower bound on bandwidth used for
%   early false returns is taken from Section 2.3 of [2]. The source code is a
%   near-exact translation of the Julia implementation in the Recognition
%   submodule of [3]. (Note that [1], [2], and [3] all use zero-based indexing
%   for bandwidth as opposed to our one-based convention.)
%
%   See also Iskcoherent.m, which uses this function to return false early in
%   some cases (minimum bandwidth up to symmetric permutation is an upper bound
%   on factor width [4] and is easier to compute).
%
%   URL: http://www.qetlab.com/has_band_k_ordering
%
%.  References
%   [1] G. M. Del Corso and G. Manzini. Finding Exact Solutions to the Bandwidth
%       Minimization Problem. Computing 62 (1999), 189–203.
%       https://doi.org/10.1007/s006070050002.
%   [2] A. Caprara and J. Salazar-González. Laying Out Sparse Graphs with
%       Provably Minimum Bandwidth. INFORMS Journal on Computing 17 (2005), no.
%       3, 356–373. https://doi.org/10.1287/ijoc.1040.0083.
%   [3] L. M. B. Varona. Luis-Varona/MatrixBandwidth.jl: Algorithms for matrix
%       bandwidth minimization and recognition. 2018.
%       https://github.com/Luis-Varona/MatrixBandwidth.jl.
%   [4] N. Johnston, S. Moein, and S. Plosker. The factor width rank of a
%       matrix. Linear Algebra and its Applications 716 (2025), 32–59.
%       https://doi.org/10.1016/j.laa.2025.03.016.

%   requires: nothing
%   author: Luis M. B. Varona (lm.varona@outlook.com)
%   package: QETLAB
%   last updated: July 19, 2025

function [has_ordering, ordering] = has_band_k_ordering(A, k)
    A_bool = offdiag_nz_support(A);

    if ~issymmetric(A_bool)
        error('A must be structurally symmetric');
    end

    if k < 0
        error('k must be a non-negative integer');
    end

    n = size(A, 1);

    if bandwidth(A) <= k
        has_ordering = true;
        ordering = 1:n;
        return;
    end

    if band_lower_bound(A_bool) > k
        has_ordering = false;
        ordering = [];
        return;
    end

    ordering_buf = zeros(n, 1);
    adj_lists = cell(n, 1);

    for i = 1:n
        adj_lists{i} = find(A_bool(:,i));
    end

    unselected_init = 1:n;
    adj_list_init = [];
    num_placed = 0;

    ordering = add_node(unselected_init, adj_list_init, num_placed);
    has_ordering = ~isempty(ordering);

    function ordering = add_node(unselected, adj_list, num_placed)
        if num_placed == n
            ordering = ordering_buf;
            return;
        end

        for candidate = unselected
            far_nodes = ordering_buf(1:num_placed-k+1);

            if ~any(A_bool(far_nodes, candidate))
                unselected_new = setdiff(unselected, candidate);

                if ~order_is_reversed(unselected_new, num_placed, candidate)
                    adj_list_new = intersect(union(adj_list, adj_lists{candidate}), unselected_new);

                    if is_compatible(adj_list_new, num_placed)
                        ordering_buf(num_placed + 1) = candidate;
                        ordering = add_node(unselected_new, adj_list_new, num_placed + 1);

                        if ~isempty(ordering)
                            return;
                        end
                    end
                end
            end
        end

        ordering = [];
    end

    function result = order_is_reversed(unselected, candidate, num_placed)
        if num_placed == 0
            first_label = candidate;
        else
            first_label = ordering_buf(1);
        end

        if isempty(unselected)
            max_last_label = candidate;
        else
            max_last_label = max(unselected);
        end

        result = (max_last_label <= first_label);
    end

    function compat = is_compatible(adj_list, num_placed)
        l = length(adj_list);

        if l >= k
            compat = false;
            return;
        end

        placed_nodes = ordering_buf(1:num_placed);
        constraints = (1:l) + num_placed;
        latest_positions = zeros(1, l);

        for j = 1:l
            neighbor = adj_list(j);
            pos = find(A_bool(placed_nodes, neighbor), true, 'first');

            if isempty(pos)
                latest_positions(j) = inf;
            else
                latest_positions(j) = pos + k - 1;
            end
        end

        latest_positions = sort(latest_positions);
        compat = all(latest_positions >= constraints);
    end
end

function band = bandwidth(A)
    [rows, cols] = ind2sub(size(A), find(A));
    dists = abs(rows - cols);

    if isempty(dists)
        band = 0;
    else
        band = max(dists) + 1;
    end
end

function lb = band_lower_bound(A)
    if all(~A(:))
        lb = 0;
        return;
    end

    A_bool = offdiag_nz_support(A);

    if all(~A_bool(:))
        lb = 1;
        return;
    end

    n = size(A_bool, 1);
    is_complete = true;
    i = 1;

    while is_complete && i < n
        j = i + 1;

        while is_complete && j <= n
            if ~A_bool(i,j)
                is_complete = false;
            end

            j = j + 1;
        end

        i = i + 1;
    end

    if is_complete
        lb = n;
        return;
    end

    D = fw_shortest_paths(A_bool);
    alpha = 1;
    gamma = n;

    for dists = D
        finite_nonzero_dists = dists(isfinite(dists) & dists > 0);

        if ~isempty(finite_nonzero_dists)
            max_dist = max(finite_nonzero_dists);
            alpha_cand = 1;
            gamma_cand = 1;

            k_hop_neighborhood_sizes = zeros(max_dist, 1);

            for k = finite_nonzero_dists
                k_hop_neighborhood_sizes(k) = k_hop_neighborhood_sizes(k) + 1;
            end

            k_hop_neighborhood_sizes = cumsum(k_hop_neighborhood_sizes) + 1;

            for k = 1:max_dist
                num_k_hop_neighbors = k_hop_neighborhood_sizes(k);
                alpha_cand = max(alpha_cand, ceil((num_k_hop_neighbors - 1) / (2 * k)) + 1);
                gamma_cand = max(gamma_cand, ceil((num_k_hop_neighbors - 1) / k) + 1);
            end

            alpha = max(alpha, alpha_cand);
            gamma = min(gamma, gamma_cand);
        end
    end

    lb = max(alpha, gamma);
end

function D = fw_shortest_paths(A)
    n = size(A, 1);
    D = zeros(n);

    for i = 1:n-1
        for j = i+1:n
            if A(i,j)
                d = 1;
            else
                d = inf;
            end

            D(i,j) = d;
            D(j,i) = d;
        end
    end

    for k = 1:n
        for i = 1:n-1
            for j = i+1:n
                if D(i,j) > D(i,k) + D(k,j)
                    d = D(i,k) + D(k,j);
                    D(i,j) = d;
                    D(j,i) = d;
                end
            end
        end
    end
end

function A_bool = offdiag_nz_support(A)
    A_bool = A ~= 0;
    A_bool(1:size(A,1)+1:end) = false;
end