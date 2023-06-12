% StaBIL with modifications -> implementation of Buckling and Optimization 
% static analysis, buckling analysis, and optimization
% Units: m, kN

% E; modulus of elasticity
% I: second moment of area
% L: length of bar/beam
% radius of the bar/beam - r1 - diagonais e r2 banzos

L = 0.100;
h = L*sqrt(3)/2;
r1=0.001;
r2=0.001;
E = 200e6;
I1 = pi*r1^4/4;
I2 = pi*r2^4/4;
EI1 = E*I1;
EI2 = E*I2;
A1 = pi*r1^2;
A2 = pi*r2^2;
EA1 = E*A1;
EA2 = E*A2;
density = 7850; % kg/m3

% Nodes=[NodID X  Y  Z] - STABIL
Nodes=  [1     0      0  0;
         2     L/2    h  0;
         3     L      0  0;
         4     L*3/2  h  0;
         5     2*L    0  0;
         6     L*5/2  h  0;
         7     3*L    0  0;          
         8     1.6*L    2*h  0];          % reference node

% Node - MATLAB - Book
% numberNodes: number of nodes
numberNodes=7;
nodeCoordinates = Nodes( : , 2:4 );

% Check the node coordinates as follows:     
figure
plotnodes(Nodes);

% Element types -> {EltTypID EltName}
Types=  {1        'beam'};

% Sections=[SecID A      ky   kz    Ixx        Iyy       Izz   yt  yb zt zb]
Sections=  [1     A1     Inf  Inf   I1         I1        I1    r1  r1  r1 r1;
            2     A2     Inf  Inf   I2         I2        I2    r2  r2  r2 r2];

% Materials=[MatID E nu];
Materials=  [1     E 0.3];              % steel



% Elements=[EltID TypID SecID MatID n1 n2 n3]
Elements=  [1     1     1     1     1  2  8;              
            2     1     1     1     2  3  8;
            3     1     2     1     3  1  8;
            4     1     1     1     3  4  8;
            5     1     2     1     4  2  8;
            6     1     2     1     3  5  8;
            7     1     1     1     5  4  8;
            8     1     1     1     5  6  8;
            9     1     2     1     6  4  8;
           10     1     2     1     5  7  8;
           11     1     1     1     7  6  8];


% generation of coordinates and connectivities - MATLAB BOOK
numberElements=11;
elementNodes = zeros(numberElements,2);
elementNodes(:, 1) = Elements(:,5);
elementNodes(:,2) = Elements(:,6);
xx=nodeCoordinates(:,1);
yy=nodeCoordinates(:,2);


% Check node and element definitions as follows:
hold('on');
plotelem(Nodes,Elements,Types);
title('Nodes and elements');

% Degrees of freedom - STABIL
% Assemble a column matrix containing all DOFs at which stiffness is
% present in the model:
DOF=getdof(Elements,Types);

% MATLAB code Finite Elements - BOOK
% GDof: global number of degrees of freedom
GDof=3*numberNodes;

% Remove all DOFs equal to zero from the vector:
%  - 2D analysis: select only UX,UY,ROTZ
%  - clamp node 1
%  - hinge at node 5
%seldof=[0.03; 0.04; 0.05; 1.00; 5.01; 5.02];
seldof=[0.03; 0.04; 0.05; 1.00; 7.02];
DOF=removedof(DOF,seldof);
prescribedDof_Stabil=[1 2 3 20]';

% Assembly of stiffness matrix K - STABIL
K=asmkm(Nodes,Elements,Types,Sections,Materials,DOF);

% stiffness matrix - MATLAB BOOK
%[stiffness]= formStiffness2Dframe(GDof,numberElements,elementNodes,numberNodes,xx,yy,EI,EA);

% boundary conditions and solution - MATLAB BOOK
%prescribedDof=[1 8 14 15]';

% Distributed loads are specified in the global coordinate system
% DLoads=[EltID n1globalX n1globalY n1globalZ ...]
DLoads=  [5     0 -10 0 0 -10 0;
          9     0 -10 0 0 -10 0];

%P=P+elemloads(DLoads,Nodes,Elements,Types,DOF); - STABIL
P=elemloads(DLoads,Nodes,Elements,Types,DOF);

% MATLAB BOOK - Force vector
%forces_book = zeros(GDof,1);
%forces_book(9) = P(2);
%forces_book(16) = P(3);
%forces_book(11) = P(8);
%forces_book(13) = P(14);
%forces_book(20) = P(15);

% Constraint equations: Constant=Coef1*DOF1+Coef2*DOF2+ ...
% Constraints=[Constant  Coef1 DOF1  Coef2 DOF2 ...]
%Constr=       [0         1     2.01  -1    3.01;
%               0         1     2.02  -1    3.02];
 
% Add constraint equations
%[K,P]=addconstr(Constr,DOF,K,P);

% Solve K * U = P  - STABIL
U=K\P;

% solution - MATLAB BOOK
%displacements=solution(GDof,prescribedDof,stiffness,forces_book);
%displacementsX = displacements(1:numberNodes);
%displacementsY = displacements(numberNodes+1:2*numberNodes);
%displacementsROT = displacements(2*numberNodes+1:3*numberNodes);

%disp(max(displacementsX));
%disp(max(displacementsY));
%disp(max(displacementsROT));
% output displacements/reactions - MATLAB BOOK
%outputDisplacementsReactions(displacements,stiffness,...
%GDof,prescribedDof)


% Plot displacements
figure
plotdisp(Nodes,Elements,Types,DOF,U,DLoads,Sections,Materials)

% The displacements can be displayed as follows:
printdisp(Nodes,DOF,U);

% Compute element forces
Forces=elemforces(Nodes,Elements,Types,Sections,Materials,DOF,U,DLoads);

% The element forces can be displayed in a orderly table:
printforc(Elements,Forces);

% Plot element forces
figure
plotforc('norm',Nodes,Elements,Types,Forces,DLoads)
title('Normal forces')

figure
plotforc('sheary',Nodes,Elements,Types,Forces,DLoads)
title('Shear forces')

figure
plotforc('momz',Nodes,Elements,Types,Forces,DLoads)
title('Bending moments')

% Plot stresses
% Calculate maximum stress - STABIL
figure
stress1 = Stress_Max('snorm',Nodes,Elements,Types,Sections,Forces,DLoads);
title('Normal stresses due to normal forces')

figure
stress2 = Stress_Max('smomzt',Nodes,Elements,Types,Sections,Forces,DLoads);
title('Normal stresses due to bending moments around z: top')

figure
stress3 = Stress_Max('smomzb',Nodes,Elements,Types,Sections,Forces,DLoads);
title('Normal stresses due to bending moments around z: bottom')

figure
stress4 = Stress_Max('smax',Nodes,Elements,Types,Sections,Forces,DLoads);
title('Maximal normal stresses')

figure
stress5 = Stress_Max('smin',Nodes,Elements,Types,Sections,Forces,DLoads);
title('Minimal normal stresses')

stress_Max_Abs_Value = max([stress1, stress2, stress3, stress4, stress5]);
X = sprintf('O maior Stress é %d MPa',stress_Max_Abs_Value/1000);
%disp(X)

% geometric stiffness matrix - formula
[Kg]= formStiffnessGeometric2Dframe(GDof,numberElements,elementNodes,xx,yy,Forces,prescribedDof_Stabil,Sections,Elements);

%critical buckling load
[EIGVEC,lambda] = eigs(K,Kg);
lambda_pcr = sort( diag(abs(lambda)) );
sprintf('O menor autovalor de flambagem é %d ', lambda_pcr(1) )
%disp( lambda_pcr(1)  )

%% Optimization

disp('Começando com a otimização, em breve o resultado estará disponível....')

%transforming functions
% Sections (function)=[SecID A             ky   kz    Ixx                 Iyy                Izz            yt    yb    zt   zb]
Sections=  @(x1,x2)       [1     pi*x1(1)^2     Inf  Inf   pi*x1(1)^4/4    pi*x1(1)^4/4        pi*x1(1)^4/4    x1(1)    x1(1)    x1(1)   x1(1)         ;
                           2     pi*x2(1)^2     Inf  Inf   pi*x2(1)^4/4    pi*x2(1)^4/4        pi*x2(1)^4/4    x2(1)    x2(1)    x2(1)   x2(1)]; 

K = @(x1,x2) asmkm(Nodes,Elements,Types,Sections(x1,x2),Materials,DOF);

U = @(x1,x2) K(x1,x2)\P;

Forces= @(x1,x2) elemforces(Nodes,Elements,Types,Sections(x1,x2),Materials,DOF,U(x1,x2),DLoads);

%Objective Function:

fun = @(x)  abs(  -abs( Stress_Max('smin',Nodes,Elements,Types,Sections(x(1), x(2)),Forces(x(1),x(2)),DLoads)/(1000)) + 250 );

%fun = @(x) calculateMass(density, Elements, Sections(x(1),x(2)) );
%fun =@(x) abs( pi()^2* 200e6 * pi*x(1)^4/4 / (0.7*L)^2 - max( abs (elemforces(Nodes,Elements,Types,Sections(x(1)),Materials,DOF,asmkm(Nodes,Elements,Types,Sections(x(1)),Materials,DOF)\elemloads(DLoads,Nodes,Elements,Types,DOF),DLoads)), [], "all" ) );
%pcr = pi()^2* 200e6 * pi*x(1)^4/4 / (1/2*L)^2;
%A=[-1, 0;0, 1];
%b=[0.001;0.010];

A=[2 -1];
b=[0];
Aeq=[];
beq=[];

x0=[0.01,0.01];
lb = [0.001,0.001];
ub = [0.010,0.010];

x_Stress = fmincon(fun,x0,A,b,Aeq,beq,lb,ub);
sprintf('O valor otimizado para o r1 é %0.4f m e para o r2 é %0.4f m', x_Stress(1), x_Stress(2)  )

 close