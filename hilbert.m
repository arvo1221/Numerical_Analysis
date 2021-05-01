%% Gauss Elimination with Maximal Column Pivoting Method
%Hilbert Matrix
% Hx = B , x = [1,1,...,1]' H is hilbert matrix
clear all
clc
format long

n = 12;
%H = 5*rand(n);
H = hilb(n);
determinant = det(H);
B = sum(H')';
S = max(H')';
x = zeros(n,1);
P = horzcat(H,B);
maxrow = zeros(1,n);
order = (1:n)';
maxcolumn = 1;
j=1;
for j=1:n
    j;
    for i=j:n
        if maxrow(j) < abs(H(j,i))/S(i)
            maxrow(j) = abs(H(j,i))/S(i);
            maxcolumn = i;
        end
    end
    if not(maxcolumn == j)
        temp = P(maxcolumn,:);
        P(maxcolumn,:) = P(j,:);
        P(j,:) = temp;
        temp2 = order(j);
        order(j) = order(maxcolumn);
        order(maxcolumn) = temp2;
    end

    for i=1:n
        for k=j+1:n
            P(k,:) = P(k,:) - (P(k,j)/P(j,j))*P(j,:);
        end
    end
    reducedP = P;
    reducedP(:,j) = [];
    reducedP(:,n-j+1) = [];
    S = max(reducedP')';
end
a = 0;
for j=n:-1:1
    if not(j==n)
       for k=j+1:n    
           a = a + x(k)*P(j,k);
       end
    end
    x(j) = (P(j,n+1) - a)/P(j,j);
    a = 0;
end
for i=1:n
    if i==1
        fprintf('x(1) is %.10f\n' ,x(i))
    elseif i==2
        fprintf('x(2) is %.10f\n' ,x(i))
    else
        fprintf('x(%d) is %.10f\n' ,i,x(i))
    end
end

%reducedP = P;
%reducedP(:,n+1) = [];
%reducedP*x - P(:,n+1);

