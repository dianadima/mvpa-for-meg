function D = crossvalED (D_train, D_test)
%cross-validated squared Euclidean distance
%input: stimulus by features matrix (computes pairwise distance between row vectors)

[n,p] = size(D_train);

k = 1;
for i = 1:n-1
    
    dsq = zeros(n-i,1);
    
    for q = 1:p
        
        dtrain = D_train(i,q) - D_train((i+1):n,q);
        dtest = D_test(i,q) - D_test((i+1):n,q);        
        dsq = dsq + dtrain.*dtest; %sum of squares replaced by sum of dot products
        
    end;
    
    D(k:(k+n-i-1)) = dsq;
    % D(k:(k+n-i-1)) = real(sqrt(dsq)); %as CV ED can be negative, the sqrt can lead to imaginary numbers
    k = k + (n-i);

end;