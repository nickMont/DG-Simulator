function nashEqbm = LH2(varargin)
    % function nashEqbm = LEMKEHOWSON(varargin)
    %
    % This function computes a sample mixed strategy Nash equilibrium in a
    % bimatrix game.  This function implements the Lemke-Howson complementary
    % pivoting algorithm for solving Bimatrix Games, a variant of the Lemke 
    % algorithm for linear complementarity problems (LCPs).
    %
    % Syntax:
    %   nashEqbm = LEMKEHOWSON(A, B)
    %   nashEqbm = LEMKEHOWSON(A, B, k0)
    %   nashEqbm = LEMKEHOWSON(A, B, k0, maxPivots)
    %
    % Parameters:
    %   A an m*n payoff matrix for the row player
    %   B an m*n payoff matrix for the column player
    %   k0 an initial pivot in the set {1,...,m+n}
    %       (optional: default = 1)
    %   maxPivots the maximum number of pivoting steps before termination
    %       (optional: default = 500000);
    %
    % Return:
    %   nashEqbm a 2x1 cell array where nashEqbm{1} and nashEqbm{2} are mixed
    %   strategies for the row and column player, respectively.
    %
    % Author:
    %	Richard M. Katzwer
    %   9/18/2013
    %   princeton.edu/~rkatzwer
    %
    % Updated:
    %   11/13/2013
    % 
    % References:
    %    C. E. Lemke and J. T. Howson, Jr.
    %      "Equilibrium Points of Bimatrix Games"
    %      Journal of the Society for Industrial and Applied Mathematics.
    %      Vol. 12, No. 2 (Jun., 1964), pp. 413-423
    %    Lloyd S. Shapley. "A note on the Lemke-Howson algorithm".
    %      Pivoting and Extension: Mathematical Programming Studies
    %      Volume 1, 1974, pp 175-189
    %    Bruno Codenotti, Stefano De Rossi, Marino Pagan.
    %      "An experimental analysis of Lemke-Howson algorithm."
    
    % Check inputs
    
    if nargin < 2 || nargin > 4
        error('This function takes between two and four arguments');
    end
    
    A = varargin{1};
    B = varargin{2};
    
    if any(size(A) ~= size(B))
		error('Matrices must have same dimension');
    end
    
    [m,n] = size(A);
    size_ = [m,n];
    
    if nargin > 2
        k0 = varargin{3};
        if k0 < 1 || k0 > m+n
            error(['Initial pivot must be in {1,...,' num2str(n+m) '}']);
        end
    else
        k0 = 1;
    end
    
    if nargin == 4
        maxPivots = varargin{4};
        if maxPivots < 1
            error('Maximum pivots parameter must be a positive integer!');
        end
    else
	    maxPivots = 500000;
    end
    
    % Scale payoffs to be strictly positive
    minVal = min( min(min(A)), min(min(B)) );
    if minVal <= 0
        A = A + ones(size(A))*(1-minVal);
        B = B + ones(size(A))*(1-minVal);
    end
    
    % Build Tableaus
    Tab = cell(2,1);
    Tab{1} = [B',     eye(n), ones(n,1)];
    Tab{2} = [eye(m), A,      ones(m,1)];
    
    % Declare row labels
    rowLabels = cell(2,1);
    rowLabels{1} = m+1:m+n;
    rowLabels{2} = 1:m;
    
    % Do complementary pivoting
    k = k0;
    if k0 <= m
        player = 1;
    else
        player = 2;
    end
    
    % Pivoting loop
    numPiv = 0;
    while numPiv < maxPivots
        
        numPiv = numPiv+1;
        
        % Use correct Tableau
        LP = Tab{player};
        [m_, ~] = size(LP);
        
        % Find pivot row (variable exiting)
        max_ = 0;
        ind = -1;
        for i = 1:m_
            t = LP(i,k) / LP(i, m+n+1);
            if t > max_
                ind = i;
                max_ = t;
            end
        end
        
        if max_ > 0
            Tab{player} = pivot(LP, ind, k);
        else
            break;
        end
        
        % swap labels, set entering variable
        temp = rowLabels{player}(ind);
        rowLabels{player}(ind) = k;
        k = temp;
        
        % If the entering variable is the same
        % as the starting pivot, break
        if k == k0
            break;
        end
        
		% update the tableau index
        if player == 1
            player = 2;
        else
            player = 1;
        end
        
    end
    
    if numPiv == maxPivots
        error(['Maximum pivot steps (' num2str(maxPivots) ') reached!']);
    end
    
    % Extract the Nash equilibrium
    nashEqbm = cell(2,1);
    
    for player = 1:2
        
        x = zeros(size_(player), 1);
        rows = rowLabels{player};
        LP = Tab{player};
        
        for i = 1:length(rows)
            if player == 1 && rows(i) <= size_(1)
                x(rows(i)) = LP(i,m+n+1) / LP(i,rows(i));
            elseif player == 2 && rows(i) > size_(1);
                x(rows(i)-size_(1)) = LP(i,m+n+1) / LP(i,rows(i));
            end
        end
        
        nashEqbm{player} = x/sum(x);
        
    end
    
end

function B = pivot(A,r,s)
% Pivots the tableau on the given row and column
    
    [m,~] = size(A);
    B = A;
    
    for i = 1 : m
        if i == r
            continue;
        else
            B(i,:) = A(i,:) - A(i,s) / A(r,s) * A(r,:);
        end
    end
    
end