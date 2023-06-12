function Q = q_dkt(b,c,det)

%Q_DKT   Q matrix for a DKT element.
%
%   Q = Q_DKT(b,c,det) returns the Q matrix of the DKT element,
%   which is used to compute the curvatures of the plate element
%
%   b       Geometrical property of the triangle (see ke_dkt) (3 * 1)
%   c       Geometrical property of the triangle              (3 * 1)
%   det     Determinant of the parametric transformation (=2*area of
%           triangle)
%
%   See also KE_DKT, KELCS_SHELL4, SE_SHELL4.
   


CC = [-b(2)/det/4 -b(3)/det/4; -c(2)/det/4 -c(3)/det/4];

p1 = 6*c(1)/(b(1)^2+c(1)^2);
p2 = 6*c(2)/(b(2)^2+c(2)^2);
p3 = 6*c(3)/(b(3)^2+c(3)^2);
q1 = 3*b(1)*c(1)/(b(1)^2+c(1)^2);
q2 = 3*b(2)*c(2)/(b(2)^2+c(2)^2);
q3 = 3*b(3)*c(3)/(b(3)^2+c(3)^2);
r1 = 3*b(1)^2/(b(1)^2+c(1)^2);
r2 = 3*b(2)^2/(b(2)^2+c(2)^2);
r3 = 3*b(3)^2/(b(3)^2+c(3)^2);
s1 = 3*c(1)^2/(b(1)^2+c(1)^2);
s2 = 3*c(2)^2/(b(2)^2+c(2)^2);
s3 = 3*c(3)^2/(b(3)^2+c(3)^2);
t1 = 6*b(1)/(b(1)^2+c(1)^2);
t2 = 6*b(2)/(b(2)^2+c(2)^2);
t3 = 6*b(3)/(b(3)^2+c(3)^2);
ro1 = (b(1)*b(2)+c(1)*c(2))/(b(1)*b(1)+c(1)*c(1));
ro2 = (b(2)*b(3)+c(2)*c(3))/(b(2)*b(2)+c(2)*c(2));
ro3 = (b(3)*b(1)+c(3)*c(1))/(b(3)*b(3)+c(3)*c(3));

AA =zeros(2,9);
for i=1:2
    e1=CC(i,1);
    e2=CC(i,2);
    AA(i,1)=e1*(30*ro2-114)-e2*(30*ro3+144);
    AA(i,4)=e1*(30*ro1+144)+e2*30*(ro1+ro3+1);
    AA(i,7)=-e1*30*(ro1+ro2+1)+e2*(114-30*ro1);
    AA(i,2)=e1*(c(2)*(14-15*ro2)-7*c(3))+e2*(7*c(2)-c(3)*(29+15*ro3));
    AA(i,3)=-e1*(b(2)*(14-15*ro2)-7*b(3))-e2*(7*b(2)-b(3)*(29+15*ro3));
    AA(i,5)= e1*(c(1)*(29+15*ro1)-7*c(3))+e2*(c(1)*(22+15*ro1)-c(3)*(15*ro3-7));
    AA(i,6)=-e1*(b(1)*(29+15*ro1)-7*b(3))-e2*(b(1)*(22+15*ro1)-b(3)*(15*ro3-7));
    AA(i,8)=e1*(c(1)*(15*ro1-7)-c(2)*(15*ro2+22))+e2*(7*c(2)+c(1)*(15*ro1-14));
    AA(i,9)=-e1*(b(1)*(15*ro1-7)-b(2)*(15*ro2+22))-e2*(7*b(2)+b(1)*(15*ro1-14));
end

G = zeros(10,9);

G(1,1)=3*p2+4*p3+AA(1,1);
G(1,2)=16-3*s2+4*s3+AA(1,2);
G(1,3)=3*q2-4*q3+AA(1,3);
G(1,4)=-3*p1-4*p3+AA(1,4);
G(1,5)=2-3*s1+4*s3+AA(1,5);
G(1,6)=3*q1-4*q3+AA(1,6);
G(1,7)=3*(p1-p2)+AA(1,7);
G(1,8)=9-3*(s1+s2)+AA(1,8);
G(1,9)=3*(q1+q2)+AA(1,9);
G(2,1)=-4*p2-3*p3+AA(1,1);
G(2,2)=16+4*s2-3*s3+AA(1,2);
G(2,3)=3*q3-4*q2+AA(1,3);
G(2,4)=3*(p3-p1)+AA(1,4);
G(2,5)=9-3*(s1+s3)+AA(1,5);
G(2,6)=3*(q1+q3)+AA(1,6);
G(2,7)=3*p1+4*p2+AA(1,7);
G(2,8)=2-3*s1+4*s2+AA(1,8);
G(2,9)=3*q1-4*q2+AA(1,9);
G(4,1)=4*(p2-p3)-AA(1,1);
G(4,2)=-9-4*(s2+s3)-AA(1,2);
G(4,3)=4*(q2+q3)-AA(1,3);
G(4,4)=10*p1+4*p3-AA(1,4);
G(4,5)=10*s1-9-4*s3-AA(1,5);
G(4,6)=4*q3-10*q1-AA(1,6);
G(4,7)=-10*p1-4*p2-AA(1,7);
G(4,8)=10*s1-9-4*s2-AA(1,8);
G(4,9)=4*q2-10*q1-AA(1,9);
G(6,1)=-3*t2-4*t3+AA(2,1);
G(6,2)=3*q2-4*q3+AA(2,2);
G(6,3)=16-3*r2+4*r3+AA(2,3);
G(6,4)=3*t1+4*t3+AA(2,4);
G(6,5)=3*q1-4*q3+AA(2,5);
G(6,6)=2+4*r3-3*r1+AA(2,6);
G(6,7)=3*(t2-t1)+AA(2,7);
G(6,8)=3*(q1+q2)+AA(2,8);
G(6,9)=9-3*(r1+r2)+AA(2,9);
G(7,1)=4*t2+3*t3+AA(2,1);
G(7,2)=3*q3-4*q2+AA(2,2);
G(7,3)=16+4*r2-3*r3+AA(2,3);
G(7,4)=3*(t1-t3)+AA(2,4);
G(7,5)=3*(q1+q3)+AA(2,5);
G(7,6)=9-3*(r1+r3)+AA(2,6);
G(7,7)=-3*t1-4*t2+AA(2,7);
G(7,8)=3*q1-4*q2+AA(2,8);
G(7,9)=2+4*r2-3*r1+AA(2,9);
G(9,1)=4*(t3-t2)-AA(2,1);
G(9,2)=4*(q2+q3)-AA(2,2);
G(9,3)=-9-4*(r2+r3)-AA(2,3);
G(9,4)=-10*t1-4*t3-AA(2,4);
G(9,5)=-10*q1+4*q3-AA(2,5);
G(9,6)=10*r1-9-4*r3-AA(2,6);
G(9,7)=10*t1+4*t2-AA(2,7);
G(9,8)=-10*q1+4*q2-AA(2,8);
G(9,9)=10*r1-9-4*r2-AA(2,9);

G(3,:)=-G(1,:);
G(5,:)=-G(2,:);
G(8,:)=-G(6,:);
G(10,:)=-G(7,:);

G(3,2)=G(3,2)+7;
G(3,5)=G(3,5)-7;
G(5,2)=G(5,2)+7;
G(5,8)=G(5,8)-7;
G(8,3)=G(8,3)+7;
G(8,6)=G(8,6)-7;
G(10,3)=G(10,3)+7;
G(10,9)=G(10,9)-7;

Q = zeros(9,9);

Q(1,:)=(b(2)*G(1,:)+b(3)*G(2,:))/7;
Q(2,:)=(2*b(2)*G(3,:)+b(3)*G(4,:))/7;
Q(3,:)=(b(2)*G(4,:)+2*b(3)*G(5,:))/7;
Q(4,:)=(c(2)*G(6,:)+c(3)*G(7,:))/7;
Q(5,:)=(2*c(2)*G(8,:)+c(3)*G(9,:))/7;
Q(6,:)=(c(2)*G(9,:)+2*c(3)*G(10,:))/7;
Q(7,:)=(c(2)*G(1,:)+c(3)*G(2,:)+b(2)*G(6,:)+b(3)*G(7,:))/7;
Q(8,:)=(2*c(2)*G(3,:)+c(3)*G(4,:)+2*b(2)*G(8,:)+b(3)*G(9,:))/7;
Q(9,:)=(c(2)*G(4,:)+2*c(3)*G(5,:)+b(2)*G(9,:)+2*b(3)*G(10,:))/7;

% formules omzetten voor theta ipv beta
% t = blkdiag(1,[0 -1;1 0]);
% T = blkdiag(t,t,t);
% Q = Q*T;
tmp = Q(:,[2 5 8]);
Q(:,[2 5 8]) = Q(:,[3 6 9]);
Q(:,[3 6 9]) = -tmp;