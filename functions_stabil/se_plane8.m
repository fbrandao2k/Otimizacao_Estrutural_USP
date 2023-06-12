function [SeGCS,SeLCS,vLCS] = se_plane8(Node,Section,Material,UeGCS,Options,gcs)
%SE_PLANE8   Compute the element stresses for a plane8 element.
%
%   [SeGCS,SeLCS,vLCS] = se_plane8(Node,Section,Material,UeGCS,Options,GCS)
%   [SeGCS,SeLCS]      = se_plane8(Node,Section,Material,UeGCS,Options,GCS)
%    SeGCS             = se_plane8(Node,Section,Material,UeGCS,Options,GCS)
%   computes the element stresses in the global and the
%   local coordinate system for the shell8 element.
%
%   Node       Node definitions           [x y z] (8 * 3)
%              Nodes should have the following order:
%              4----7----3
%              |         |
%              8         6
%              |         |
%              1----5----2
%   Section    Section definition         [h]  (only used in plane stress)
%   Material   Material definition        [E nu rho]
%   UeGCS      Displacements (8 * nTimeSteps)
%   Options    Element options            {Option1 Option2 ...}
%   GCS        Global coordinate system in which stresses are returned
%              'cart'|'cyl'
%   SeGCS      Element stresses in GCS in corner nodes IJKL
%              48 = 6 stress comp. * 8 nodes (24 * nTimeSteps)
%                                        [sxx syy szz sxy syz sxz]
%   SeLCS      Element stresses in LCS in corner nodes IJKL
%              48 = 6 stress comp. * 8 nodes (24 * nTimeSteps)
%                                        [sxx syy szz sxy syz sxz]
%   vLCS       Unit vectors of LCS (1 * 9)
%
%   See also ELEMSTRESS, SE_SHELL8.

% Stijn Fran�ois
% 2016

% Options
if nargin<5, Options=[]; end
if ~isfield(Options,'problem'), Options.problem='2dstress'; end

% Check nodes
Node=Node(1:8,1:3);

SeLCS = selcs_plane8(Node,Section,Material,UeGCS,Options);

if ~isempty(gcs)
  switch lower(gcs)
  case  'cart'
    SeGCS=SeLCS;
  case  'cyl'
    a1 = [Node(1:8,1) Node(1:8,2) zeros(4,1)];
    a2 = [-Node(1:8,2) Node(1:8,1) zeros(4,1)];
  
     for iNode=1:8
     A= [a1(iNode,:)/norm(a1(iNode,:));a2(iNode,:)/norm(a2(iNode,:));0 0 1];
     theta = vtrans_solid(A);
        for z = 1:3
        SeGCS((1:6)+6*(iNode-1)+24*(z-1),:)= theta*SeGCS((1:6)+6*(iNode-1)+24*(z-1),:);
        end
     end
  case 'sph'
    error('Sperhical coordinate system not defined for 2D problem.');        
  end
end

vLCS=[];
end
