function [Ke,Me,dKedx]=ke_plane15(Node,Section,Material,Options,dNodedx,dSectiondx)

%KE_PLANE15   plane element stiffness and mass matrix in global coordinate system.
%
%   [Ke,Me] = ke_plane15(Node,Section,Material,Options) returns the element
%   stiffness and mass matrix in the global coordinate system for a 15-node
%   plane triangular element. Plane15 only operates in the 2D xy-plane so that
%   z-coordinates should be equal to zero.
%
%   Node       Node definitions           [x y z] (15 * 3)
%              Nodes should have the following order:
%
%                  3
%                  | \
%                 10   9
%                  |     \
%                 11  15   8
%                  |         \
%                 12  13   14  7
%                  |             \
%                  1---4---5---6---2
%
%   Section    Section definition         [h] (only used in plane stress)
%   Material   Material definition        [E nu rho]
%   Options    Struct containing optional parameters. Fields:
%      .problem Plane stress or plane strain
%               {'2dstress' (default) | '2dstrain'}
%      .nXi     Number of Gauss integration points
%               {1 | 3 | 4 | 6 | 7  | 12 | 13 (default) | 16 | 19 | 28 | 33 | 37}
%   Ke         Element stiffness matrix (30 * 30)
%   Me         Element mass matrix (30 * 30)
%
%   See also KE_BEAM, ASMKM, KE_TRUSS.

% Stijn Fran�ois
% 2017

% Options
if nargin<4, Options=[]; end
if ~isfield(Options,'problem'), Options.problem='2dstress'; end
if ~isfield(Options,'nXi'),     Options.nXi=13; end

if nargin<5, dNodedx = []; end
if nargin<6, dSectiondx = []; end
if (nnz(dNodedx)+nnz(dSectiondx)) > 0
    error('Sensitivities have not been implemented yet.')
end
nVar = 0;
if nargout>2 
    if ~isempty(dNodedx) || ~isempty(dSectiondx)
        nVar = max(size(dNodedx,3),size(dSectiondx,3));
    end
    dKedx = cell(nVar,1);
end

if nargout>1, Me = []; end
if (nargout==2), rho=Material(3); end

% Constitutive matrix
C=cmat_isotropic(Options.problem,Section,Material);

% Gauss integration points
[xi,H]=gaussqtri(Options.nXi);

Ke=zeros(30);
if (nargout==2), Me=zeros(30); end

for iXi=1:Options.nXi
    [Ni,dNi_dxi,dNi_deta]=sh_t15(xi(iXi,1),xi(iXi,2));

    J=[dNi_dxi dNi_deta].'*Node(1:15,1:2);
    detJ=det(J);

    dN=J\[dNi_dxi dNi_deta].';

    B=zeros(3,30);
    B(1,1:2:29)=dN(1,:);
    B(2,2:2:30)=dN(2,:);
    B(3,1:2:29)=dN(2,:);
    B(3,2:2:30)=dN(1,:);

    Ke=Ke+B.'*C*B*H(iXi)*detJ;

    if (nargout==2)
    	N=zeros(2,30);
        N(1,1:2:29)=Ni;
        N(2,2:2:30)=Ni;
        Me=Me+N.'*rho*N*H(iXi)*detJ;
    end
end
end
