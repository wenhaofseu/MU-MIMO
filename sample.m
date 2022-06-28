function location = sample(radius,center,K)
    r = sqrt(radius^2*rand(K,1));
    theta = 2*pi*rand(K,1);
    location = repmat(center,[K,1]) +[r.*cos(theta),r.*sin(theta),zeros(K,1)];
end

