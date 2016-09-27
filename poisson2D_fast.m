% FEM Study Code - Poisson's Equation (-U_xx - U_yy = f) on a Unit Square
% Author: Casey Tyler Icenhour
% Organization: North Carolina State Univ./Oak Ridge National Laboratory
% September 2016

% Inspired by 2D FEM Poisson code by Gilbert Strang at MIT in his text:
% "Computational Science and Engineering" 2nd Edition, Wellesley-Cambridge 
% Press, 2008.
% See http://math.mit.edu/~gs/cse/ for more info
% Relevant MATLAB Code: http://math.mit.edu/~gs/cse/codes/femcode.m

clear all
close all

% Generate 2D mesh using DistMesh (see http://persson.berkeley.edu/distmesh/
% for more info)
% KEY NOTES FOR DISTMESH
% p = node positions -- N-by-2 array contains x,y coordinates for each of
%     the N nodes
% t = triangle indices -- row associated with each triangle has 3 integer
%     entries to specify node numbers in the triangle 

% Generate 2D mesh on unit square
len = 1;
wid = 1;
nom = wid/5;   % nominal edge width of triangles (DistMesh uses for 
                % generation, is more than likely smaller than this)

% Specify distance function for geometry (rectangular waveguide)
fd=@(p) drectangle(p,0,wid,0,len);

% Note: @huniform notes uniform edge length for triangles
[p,t]=distmesh2d(fd,@huniform,nom,[0,0;wid,len],[0,0;wid,0;0,len;...
    wid,len]);

% Generate list of boundary edges for BC using DistMesh functions

e = boundedges(p,t);    % N-by-2 array - each row is an edge defined by 
                        %                two nodes
b = unique(e);          % column vector of edge nodes

% number of nodes, number of triangles
N = size(p,1);
T = size(t,1);

% Initialize K matrix and F vector in system KU=F, where U is the solution
K = sparse(N,N); % Sparse format
F = zeros(N,1);

% Integrating the weak form via an element-oriented assembly method (i.e.
% over each triangle element at a time in this case)
% NOTE: Using "hat" nodal test functions and one-point Gauss quadrature for
% integrations
for i=1:T 
    nodes = t(i,:); % pull relevant nodes of triangle
    
    % REALLY GREAT TRICK FOR FINDING AREA OF A TRIANGLE USING
    % DETERMINANTS (aka matrices are awesome)
    Pe = [ones(3,1),p(nodes,:)];    % each row is: 1  x_coord  y_coord
    Area = abs(det(Pe))/2;
    
    % Since we're using linear test/basis functions (phi), form of 
    % functions inside triangle is a + bx + cy = phi. For a given phi tied 
    % to a particular node, we want it to be 1 at the node and 0 at all
    % surrounding nodes in the triangle. We can use Pe to get this, since
    % the relevant system of equations is Pe * coeff = I. The columns of
    % the matrix coeff are the coefficients of the function form of phi at
    % the relevant node. For example, a1, b1, and c1 is column 1 of coeff
    % referring to the function phi_1 at node 1 of the current triangle
    % element. So, inverting Pe will give us coeff.
    
    coeff = inv(Pe); 
    
    % Using coeff, getting the derivatives of phi in x and y is 
    % easy - we only need to take the bottom two rows of coeff
    
    grad = coeff(2:3,:);
    
    % We need to form the element stiffness matrix (see Gockenbach's 2006 
    % FEM book pg. 130 for details). This is the component of matrix K
    % relevant to the triangle under inspection. 
    %
    % One point Gauss quadrature for triangles says that each integral
    % can be calculated by using: 
    % [Trangle Area]*[Grad of test function]*[Grad of basis function]
    % Where the test function and basis function are the parts of the same
    % basis in the Galerkin formulation
    
    Ke = Area*grad'*grad;
    
    % Since f(x,y) = 1, the integral over the test function is just the
    % area of the triangular pyramid of height 1
    Fe = Area/3;
    
    % Assemble K and F
    
    K(nodes,nodes) = K(nodes,nodes) + Ke;
    F(nodes) = F(nodes) + Fe;
end

% Dirichlet Boundary Conditions: U(x,y) = 0 at all boundary nodes (in b)
% Set all boundary rows and columns of K and F to be zero
K(b,:) = 0;
% K(:,b) = 0;
F(b) = 0;

% Insert identity matrix into boundary submatrix of K
% WHY?
%   NOTE: Without this line, the resulting matrix is singular, according
%         to MATLAB - adjusting it to another number (0.5, 2, and 5 tested) 
%         leads to a singular matrix as well 
%   Considering the structure of K (row/column positions corresponding
%   to basis and test functions in the integrals), this seems to just allow
%   for unweighted inclusion of the U(x,y) boundary terms in the
%   discretized solution for a particular position at the boundary. (i.e.
%   for boundary node 1, having K(1,1) = 1 allows for U(1,1) to be
%   represented in the resulting first equation in the sytem (here, it
%   seems, that is U(1,1) = 0). Otherwise some of the rows of K are
%   trivial (i.e. 0 = 0).

K(b,b) = speye(length(b),length(b));

% Solve system of equations 
U = K\F;

% Plot solution
trisurf(t,p(:,1),p(:,2),0*p(:,1),U,'edgecolor','k','facecolor','interp');
view(2),axis equal,colorbar
title('FEM solution')

% Generate and plot analytic solution on triangular mesh
U_analytic = analyticFXN_poisson2D(p,50);

figure
trisurf(t,p(:,1),p(:,2),0*p(:,1),U_analytic,'edgecolor','k','facecolor','interp')
view(2),axis equal,colorbar
title('Analytic solution to 50 terms in infinite series');

two_norm_error = norm(U_analytic-U,2)

% figure
% trisurf(t,p(:,1),p(:,2),0*p(:,1),U_analytic-U,'edgecolor','k','facecolor','interp')
% view(2),axis equal,colorbar
% title('Absolute error')