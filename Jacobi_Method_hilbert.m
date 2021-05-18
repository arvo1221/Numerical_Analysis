clear all
clc

MaxIter = 10;
TOL = 1e-04;

n = 4;
H = hilb(n);
B = sum(H')';

x0 = zeros(n,1);
x_old = x0;
At = zeros(n,n); bt = zeros(n,1);
for i=1:n
   At(i,:) = H(i,:)/H(i,i);
   At(i,i) = 0;
   bt(i) = (B(i)/H(i,i));
end

for repeat=1:MaxIter
    x = -At*x_old + bt;
    if norm(x-x_old,2)<TOL
        break;
    end
    size = norm(x-x_old,2);
    x_old = x;
end

for i=1:n
    
    if i==1
        fprintf('x(1) is %.10f\n' ,x_old(i))
    elseif i==2
        fprintf('x(2) is %.10f\n' ,x_old(i))
    else
        fprintf('x(%d) is %.10f\n' ,i,x_old(i))
    end
end
fprintf('error = %.2f\n',size)
fprintf('repeat %d\n' ,repeat)