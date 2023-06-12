function  [Kg]=...
    formStiffnessGeometric2Dframe(GDof,numberElements,...
    elementNodes,xx,yy,Forces, removeDOFs,Sections, Elements)

Kg=zeros(GDof); 
% computation of the geometric stiffness matrix
for e=1:numberElements
  %get the Area of the element
  A = Sections( Elements(e,3), 2);
  % elementDof: element degrees of freedom (Dof)
  indice=elementNodes(e,:)   ;       
  elementDof=[ indice*3-2 indice*3-1 indice*3] ;
  %nn=length(indice);  
  xa=xx(indice(2))-xx(indice(1));
  ya=yy(indice(2))-yy(indice(1));  
  length_element=sqrt(xa*xa+ya*ya);
  cosa=xa/length_element;  
  sena=ya/length_element;
  ll=length_element;
  
  L= [cosa*eye(2) sena*eye(2) zeros(2);
     -sena*eye(2) cosa*eye(2) zeros(2);
     zeros(2,4) eye(2)];

   kg1=-Forces(e,1)/(30*ll*A)*[0 0    0        0      0       0;
                            0 36   3*ll     0      -36     3*ll;
                            0 3*ll 4*ll^2   0      -3*ll   -ll^2;
                            0 0    0        0      0       0;  
                            0 -36  -3*ll    0      36      -3*ll;
                            0 3*ll -ll^2    0      -3*ll   4*ll^2];

       Kg(elementDof,elementDof)=...
        Kg(elementDof,elementDof)+L'*kg1*L;
end

% remove Gdofs of boundary condition 
Kg_active = Kg;
j = 0;
for i = 1:length(removeDOFs)
    Kg_active(:,removeDOFs(i-j)) = [];
    Kg_active(removeDOFs(i-j), :) = [];
    j = j + 1;
end

Kg = Kg_active;



