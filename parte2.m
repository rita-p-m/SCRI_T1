format long
x=[1; 1; 1; 0; 0]
W=diag([1/0.0002; 1/0.0002; 1/0.0002; 1/0.001; 1/0.001; 1/0.001; 1/0.001; 1/0.005; 1/0.005; 1/0.00001])
z = [1.08; 1.03; 1; 11.88; 4.98; 8.3; 8; 3.2; 3; 0] %remover Q13
g = [3 -2 -1; -2 2 0; -1 0 1]
b = [-27 18 9; 18 -18 0; 9 0 -9]
y = g + 1i*b;
tol = 10; %tive de atribuir um valor senão não imprimia o que está dentro do while
a=0;
rN = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0]
while tol>10^-4
    v1=x(1)
    v2=x(2)
    v3=x(3)
    theta2=x(4)
    theta3=x(5)
    H=[1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 2*v1*g(1)+v2*(g(2)*cos(theta2)+b(2)*(-sin(theta2)))+v3*(g(3)*cos(theta3)+b(3)*(-sin(theta3))) v1*(g(2)*cos(theta2)+b(2)*(-sin(theta2))) v1*(g(3)*cos(theta3)+b(3)*sin(-theta3)) v1*v2*(g(2)*(-sin(theta2))-b(2)*cos(theta2)) v1*v3*(g(3)*(-sin(theta3))-b(3)*cos(theta3)); v3*(g(7)*cos(theta3)+b(7)*sin(theta3)) 0 v1*(g(7)*cos(theta3)+b(7)*sin(theta3))-2*g(7)*v3 0 v3*v1*(-g(7)*sin(theta3)+b(7)*cos(theta3)); v2*((g(2)/2)*cos(theta2)+(b(2)/2)*(-sin(theta2)))-g(2)*v1 v1*((g(2)/2)*cos(theta2)+(b(2)/2)*(-sin(theta2))) 0 v1*v2*((g(2)/2)*(-sin(theta2))-(b(2)/2)*cos(theta2)) 0; v2*((g(2)/2)*cos(theta2)+(b(2)/2)*(-sin(theta2)))-g(2)*v1 v1*((g(2)/2)*cos(theta2)+(b(2)/2)*(-sin(theta2))) 0 v1*v2*((g(2)/2)*(-sin(theta2))-(b(2)/2)*cos(theta2)) 0; v2*((g(2)/2)*(-sin(theta2))-(b(2)/2)*cos(theta2))+b(2)*v1 v1*((g(2)/2)*(-sin(theta2))-(b(2)/2)*cos(theta2)) 0 -(v1*v2*((g(2)/2)*cos(theta2)+(b(2)/2)*(-sin(theta2)))) 0; v2*((g(2)/2)*(-sin(theta2))-(b(2)/2)*cos(theta2))+b(2)*v1 v1*((g(2)/2)*(-sin(theta2))-(b(2)/2)*cos(theta2)) 0 -(v1*v2*((g(2)/2)*cos(theta2)+(b(2)/2)*(-sin(theta2)))) 0; v3*(g(7)*sin(theta3)-b(7)*cos(theta3)) 0 v1*(g(7)*sin(theta3)-b(7)*cos(theta3))+2*v3*b(7) 0 v3*v1*(g(7)*cos(theta3)+b(7)*sin(theta3))]
    Ht = transpose(H)
    G = Ht*W*H
    Ginv = G^-1
    h = [v1; v2; v3; v1*(v1*g(1)+v2*(g(2)*cos(theta2)+b(2)*(-sin(theta2)))+v3*(g(3)*cos(theta3)+b(3)*(-sin(theta3)))); v3*v1*(g(7)*cos(theta3)+b(7)*sin(theta3))-g(7)*v3*v3; v1*v2*((g(2)/2)*cos(theta2)+(b(2)/2)*(-sin(theta2)))-(g(2)/2)*v1*v1; v1*v2*((g(2)/2)*cos(theta2)+(b(2)/2)*(-sin(theta2)))-(g(2)/2)*v1*v1; v1*v2*((g(2)/2)*(-sin(theta2))-(b(2)/2)*cos(theta2))+v1*v1*(b(2)/2); v1*v2*((g(2)/2)*(-sin(theta2))-(b(2)/2)*cos(theta2))+v1*v1*(b(2)/2); v3*v1*(g(7)*sin(theta3)-b(7)*cos(theta3))+v3*v3*b(7)]
    r = z - h
    deltax = Ginv*Ht*W*r
    x=x+deltax
    tol = max(abs(deltax))
    a = a+1 %numero de iterações
end

zht = transpose(z-h);
J=zht*W*(z-h)

ri=z-h
R = W^-1
c=R-H*Ginv*Ht
cii=diag(c)

k=1;
while k<=10
    rN(k) = ri(k)/sqrt(cii(k))
    k=k+1;
end
