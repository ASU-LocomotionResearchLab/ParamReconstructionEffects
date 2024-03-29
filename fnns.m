function [fnn] = fnns (x, m, D, n, p) 
% FNNS  this function determines the embedding dimension for a scalar time 
%       series, using the method of false nearest neighbors similar to:
%           "Determining embedding dimension for phase-space 
%           reconstruction using a geometrical construction", 
%           M. B. Kennel, R. Brown, and H.D.I. Abarbanel,
%           Physical Review A, Vol 45, No 6, 15 March 1992, 
%           pp 3403-3411.
%       method was further improved to account for temporal correlations, 
%       autocorrelations and sparsity of data.
%
% CALL: [fnn] = fnns (x, m, D, n, p);
%
% INPUTS:
%   x: the vector containing the scalar observations
%   m: delay parameter used in embedding
%   D: maximum embedding dimension to consider
%   n: number of data points to use
%   p: (optional) length of temporal correlation to remove
%      use 0 for not removing temporal correlations, otherwise
%      default value equal to 2*m will be used.
%
% copiright (c) and written by David Chelidze 2/27/2009
 
if nargin < 4
    error('need at least four input variables')
end
if nargin < 5
    p = 2*m; % remove temporal correlations of this length
    % NOTE: if data is highly oversampled this may need to be increased
end
 
% check data and make a raw vector of scalar data
if min(size(x)) ~= 1
   error('first input should be a vector of scalar data points.')
else
    x=x(:)'; % just the way we need it
end
 
nm = length(x) - D*m; % maximum number of points possible
if n > nm
    fprintf('WARNING: Maximum number of embedded vectors exceeded\n')
    fprintf('ACTION:  continuing with maximum number possible!\n')
    n = nm
end
 
R_tol = 15; % tolerance suggested by KBA (15-20) -- can be changed!
R_a = var(x); % estimate of attractor size, as per KBA
A_tol = 4; % suggested by KBA (2-4) -- try other values!
 
indx = 1:n; % index of points in the embedding
y = []; % initialize array to hold currently embedded vectors
z = x(indx); % initialize array used for d+1 dimension coordinates
fnn = zeros(1,D); % initialize array to hold fnn list
 
for d = 1:D
    q = y; % (d-1)-dimensional embedding
    y = [q; z]; % d-dimensional embedding
    z = x(indx + m*d); % additional coordinates for (d+1)-dimensional 
    L = zeros(1,n); % initialize array for false nns
%     fprintf('Partitioning data for d = %d\n', d)
    [kd_tree, r] = kd_partition(y, 256); % partition into kd-tree
%     fprintf('Checking for false nearest neighbors\n')
    for k = 1:n % find the false nns
        s = []; % holds indices of the temporarily uncorrelated nns
        ms = 10; % try if this can be found in the first 20 nns
        while isempty(s) % search for more nns
            ms = 2*ms;
            [nns, nnd, nni] = kd_search(y(:,k), ms, kd_tree, y(:,r));
            s = find( abs(r(nni{1}) - r(nni{1}(1))) > p );
        end
        ir = r(nni{1}(s(1)));  % index of temporarily uncorrelated nn in y
        diff = abs(z(ir) - z(k)); % additional length due to d+1 dim
        if diff > nnd{1}(s(1))*R_tol || ... % check if nn moves too far
           (nnd{1}(s(1))^2 + diff^2)/R_a > A_tol % check if dim is too big
           L(k) = 1; 
        end
    end % search
    fnn(d) = sum(L)/n; % count and normalize FNNs
end % done
 
if nargout == 0 % plot the results
    plot(1:D,100*fnn,'o-') 
    xlabel('Embedding Dimension (size)')
    ylabel('FNNs (%)')
end


function [kd_tree, r] = kd_partition(y, b, c)
% KD_PARTITION  Create a kd-tree and partitioned database for efficiently 
%               finding the nearest neighbors to a point in a 
%               d-dimensional space.
%
% USE: [kd_tree, r] = kd_partition(y, b, c);
%
% INPUT:
%   y: original multivariate data (points arranged columnwise, size(y,1)=d)
%   b: maximum number of distinct points for each bin (default is 100)
%   c: minimum range of point cluster within a final leaf (default is 0)
%
% OUTPUT:
%   kd_tree structure with the following variables:
%       splitdim: dimension used in splitting the node
%       splitval: corresponding cutting point
%       first & last: indices of points in the node
%       left & right: node #s of consequent branches to the current node
%   r: sorted index of points in the original y corresponding to the leafs
%
% to find k-nearest neighbors use kd_search.m
%
% copyrighted (c) and written by David Chelidze, January 28, 2009.
 
% check the inputs
if nargin == 0
    error('Need to input at least the data to partition')
elseif nargin > 3
    error('Too many inputs')
end
 
 
% initializes default variables if needed
if nargin < 2
    b = 100;
end
if nargin < 3
    c = 0;
end
 
[d, n] = size(y); % get the dimension and the number of points in y
 
r = 1:n; % initialize original index of points in y
 
% initializes variables for the number of nodes and the last node
node = 1;
last = 1;
 
% initializes the first node's cut dimension and value in the kd_tree
kd_tree.splitdim = 0;
kd_tree.splitval = 0;
 
% initializes the bounds on the index of all points
kd_tree.first = 1;
kd_tree.last = n;
 
% initializes location of consequent branches in the kd_tree
kd_tree.left = 0;
kd_tree.right = 0;

while node <= last % do until the tree is complete
    
    % specify the index of all the points that are partitioned in this node
    segment = kd_tree.first(node):kd_tree.last(node);
    
    % determines range of data in each dimension and sorts it
    [rng, index] = sort(range(y(:,segment),2));
    
    % now determine if this segment needs splitting (cutting)
    if rng(d) > c && length(segment)>= b % then split
        yt = y(:,segment); 
        rt = r(segment);
        [sorted, sorted_index] = sort(yt(index(d),:));
        % estimate where the cut should go
        lng = size(yt,2);
        cut = (sorted(ceil(lng/2)) + sorted(floor(lng/2+1)))/2;
        L = (sorted <= cut); % points to the left of cut
        if sum(L) == lng % right node is empty
            L = (sorted < cut); % decrease points on the left
            cut = (cut + max(sorted(L)))/2; % adjust the cut
        end
        
        % adjust the order of the data in this node
        y(:,segment) = yt(:,sorted_index); 
        r(segment) = rt(sorted_index);
 
        % assign appropriate split dimension and split value
        kd_tree.splitdim(node) = index(d);
        kd_tree.splitval(node) = cut;
        
        % assign the location of consequent bins and 
        kd_tree.left(node) = last + 1;
        kd_tree.right(node) = last + 2;
        
        % specify which is the last node at this moment
        last = last + 2;
        
        % initialize next nodes cut dimension and value in the kd_tree
        % assuming they are terminal at this point
        kd_tree.splitdim = [kd_tree.splitdim 0 0];
        kd_tree.splitval = [kd_tree.splitval 0 0];
 
        % initialize the bounds on the index of the next nodes
        kd_tree.first = [kd_tree.first segment(1) segment(1)+sum(L)];
        kd_tree.last = [kd_tree.last segment(1)+sum(L)-1 segment(lng)];
 
        % initialize location of consequent branches in the kd_tree
        % assuming they are terminal at this point
        kd_tree.left = [kd_tree.left 0 0];
        kd_tree.right = [kd_tree.right 0 0];
        
    end % the splitting process
 
    % increment the node
    node = node + 1;
 
end % the partitioning

function [pqr, pqd, pqi] = kd_search(y,r,tree,yp)
% KD_SEARCH     search kd_tree for r nearest neighbors of point yq.
%               Need to partition original data using kd_partition.m
%
% USE: [pqr, pqd, pqi] = kd_search(y,r,tree,yp);
%
% INPUT:
%       y: array of query point in columnwise form
%       r: requested number of nearest neighbors to the query point yq
%       tree: kd_tree constructed using kd_partition.m
%       yp: partitioned (ordered) set of data that needs to be searched
%          (using my kd_partirion you want to input ym(:,indx), where ym 
%           is the data used for partitioning and indx sorted index of ym)
%
% OUTPUT:
%       pqr: cell array of the r nearest neighbors of y in yp 
%       pqd: cell array of the corresponding distances
%       pqi: cell array of the indices of r nearest neighbors of y in yp 
%
% copyright (c) and written by David Chelidze, February 02, 2009.
 
% check inputs
if nargin < 4
    error('Need four input variables to work')
end
 
% declare global variables for subfunctions
global yq qri qr qrd finish b_lower b_upper
 
[d, n] = size(y);
pqr = cell(n,1);
pqd = cell(n,1);
pqi = cell(n,1);
 
for k = 1:n,
    yq = y(:,k);
    qrd = Inf*ones(1,r); % initialize array for r distances
    qr = zeros(size(yq,1),r); % initialize r nearest neighbor points
    qri = zeros(1,r); % initialize index of r nearest neighbors
    finish = 0; % becomes 1 after search is complete
 
    % set up the box bounds, which start at +/- infinity (whole kd space)
    b_upper = Inf*ones(size(yq));
    b_lower = -b_upper;
 
    kdsrch(1,r,tree,yp); % start the search from the first node
    pqr{k} = qr;
    pqd{k} = sqrt(qrd);
    pqi{k} = qri;
end
 
%--------------------------------------------------------------------------            
function kdsrch(node,r,tree,yp)
% KDSRCH    actual kd search
% this drills down the tree to the end node and updates the 
% nearest neighbors list with new points
%
% INPUT: starting node number, and kd_partition data
%
% copyright (c) and written by David Chelidze, February 02, 2009.
 
global yq qri qr qrd finish b_lower b_upper
 
% first find the terminal node containing yq
if tree.left(node) ~= 0 % not a terminal node, search deeper
 
    % first determin which child node to search
    if yq(tree.splitdim(node)) <= tree.splitval(node)
        % need to search left child node
        tmp = b_upper(tree.splitdim(node));
        b_upper(tree.splitdim(node)) = tree.splitval(node);
        kdsrch(tree.left(node),r,tree,yp);
        b_upper(tree.splitdim(node)) = tmp;
    else % need to search the right child node
        tmp = b_lower(tree.splitdim(node));
        b_lower(tree.splitdim(node)) = tree.splitval(node);
        kdsrch(tree.right(node),r,tree,yp);
        b_lower(tree.splitdim(node)) = tmp;
    end % when the terminal node (leaf) containing yq is reached
    if finish % done searching
        return
    end
 
    % check if other nodes need to be searched
    if yq(tree.splitdim(node)) <= tree.splitval(node)
        tmp = b_lower(tree.splitdim(node));
        b_lower(tree.splitdim(node)) = tree.splitval(node);
        if overlap(yq, b_lower, b_upper, qrd(r)) 
            % need to search the right child node
            kdsrch(tree.right(node),r,tree,yp);
        end
        b_lower(tree.splitdim(node)) = tmp;
    else
        tmp = b_upper(tree.splitdim(node));
        b_upper(tree.splitdim(node)) = tree.splitval(node);
        if overlap(yq, b_lower, b_upper, qrd(r)) 
            % need to search the left child node
            kdsrch(tree.left(node),r,tree,yp);
        end
        b_upper(tree.splitdim(node)) = tmp;
    end % when all the other nodes are searched
    if finish % done searching
        return
    end
 
else % this is a terminal node: update the nearest neighbors
 
    yt = yp(:,tree.first(node):tree.last(node)); % all points in node
    dstnc = zeros(1,size(yt,2));
    for k = 1:size(yt,1)
        dstnc = dstnc + (yt(k,:)-yq(k)).^2; 
    end
    %dstnc = sum((yt-yq*ones(1,size(yt,2))).^2,1); % distances squared
    qrd = [qrd dstnc]; % current list of distances squared
    qr = [qr yt]; % current list of nearest neighbors
    qri = [qri tree.first(node):tree.last(node)];
    [qrd, indx] = sort(qrd); % sorted distances and their index
    qr = qr(:,indx); % sorted list of nearest neighbors
    qri = qri(indx); % sorted list of indexes
    if size(qr,2) > r % truncate to the first r points
        qrd = qrd(1:r);
        qr = qr(:,1:r);
        qri = qri(1:r);
    end
    % be done if all points are with this box on the first run
    if within(yq, b_lower, b_upper, qrd(r));
        finish = 1;
    end % otherwise (during backtracking) WITHIN will always return 0.
    
end % if
 
%--------------------------------------------------------------------------            
function flag = within(yq, b_lower, b_upper, ball)
% WITHIN    check if additional nodes need to be searched (i.e. if the ball
%  centered at yq and containing all current nearest neighbors overlaps the
%  boundary of the leaf box containing yq)
%
% INPUT:
%   yq: query point
%   b_lower, b_upper: lower and upper bounds on the leaf box
%   ball: square of the radius of the ball centered at yq and containing
%         all current r nearest neighbors
% OUTPUT:
%   flag: 1 if ball does not intersect the box, 0 if it does
%
% Modified by David Chelidze on 02/03/2009.
 
if ball <= min([abs(yq-b_lower)', abs(yq-b_upper)'])^2
    % ball containing all the current nn is inside the leaf box (finish)
    flag = 1;
else % ball overlaps other leaf boxes (continue recursive search)
    flag = 0; 
end
 
%--------------------------------------------------------------------------            
function flag = overlap(yq, b_lower, b_upper, ball)
% OVERLAP   check if the current box overlaps with the ball centered at yq
%   and containing all current r nearest neighbors�.
%
% INPUT:
%   yq: query point
%   b_lower, b_upper: lower and upper bounds on the current box
%   ball: square of the radius of the ball centered at yq and containing
%         all current r nearest neighbors
% OUTPUT:
%   flag: 0 if ball does not overlap the box, 1 if it does
%
% Modified by David Chelidze on 02/03/2009.
 
il = find(yq < b_lower); % index of yq coordinates that are lower the box 
iu = find(yq > b_upper); % index of yq coordinates that are upper the box
% distance squared from yq to the edge of the box
dst = sum((yq(il)-b_lower(il)).^2,1)+sum((yq(iu)-b_upper(iu)).^2,1);
if dst >= ball % there is no overlap (finish this search)    
    flag = 0;
else % there is overlap and the box needs to be searched for nn
    flag = 1;
end