clear all
clc
format long

MaxIter = 1000000;
TOL = 0.99899e-04;
n = 15;
H = hilb(n);
B = sum(H')';

x0 = zeros(n,1);
%{
when omega = 1 , Gauss-Seidel Method
when 1<omegaa<2 , successive over-relaxation method
%}
omega = 1.5;
x_old = x0;
x = x0;

repeat = 1;
for repeat=1:MaxIter
    x(1) = x(1) - omega*(H(1,:)*x - B(1))/H(1,1);
    for i = 2:n-1
       x(i) = x(i) - omega/H(i,i)*(H(i,1:i-1)*x(1:i-1)+H(i,i:n)*x(i:n)-B(i)); 
    end
    x(n) = x(n) - omega*(H(n,:)*x - B(n))/H(n,n);
    if norm(x-x_old,2) < TOL
        break;
    end
    size = norm(x-x_old,2);
    x_old = x;
end

for i=1:n
    if i==1
        fprintf('x(1) is %.2f   x(1) abs error is %.2f \n' ,x(i),x(i)-1)
    elseif i==2
        fprintf('x(2) is %.2f   x(2) abs error is %.2f \n' ,x(i),x(i)-1)
    else
        fprintf('x(%d) is %.2f   x(%d) abs error is %.2f \n' ,i,x(i),i,x(i)-1)
    end
end
fprintf('error = %f\n',size)
fprintf('repeat %d\n ' ,repeat)