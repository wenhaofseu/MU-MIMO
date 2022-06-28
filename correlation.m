function matrix = correlation(rho,N)
    matrix = zeros(N);
    for i=1:N
        for j=1:N
            matrix(i,j) = rho^abs(i-j);
        end
    end
end

