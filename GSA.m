% this function is used for part 1 only - Gram Shmidt Algorithm

function [x1, x2] = GSA(s, s_ij, phi_j, fs, M)
    E_1 = sum(s(1, :).^2)*(1/fs); % Es1 
    if E_1 == 0
        phi_j(1, :) = 0; % phi_1(t) = s1(t)/Es1
    else
        phi_j(1, :) = s(1, :) / sqrt(E_1);
    end
    s_ij(1, 1) = sqrt(E_1);
    
    for i=2:M
        phi_s = [zeros(1, length(phi_j(1, :)))];
        for j=1:i-1
            s_ij(i, j) = sum(s(i, :) .* phi_j(j, :)) * (1/fs);
            phi_s = phi_s + s_ij(i, j)*phi_j(j, :);
        end
        g_i = s(i, :) - phi_s;
        if s(i, :) == phi_s
            phi_j(i, :) = 0;
        else
            Eg = sum(g_i.^2)*(1/fs);
            s_ij(i, i) = sum(s(i, :).*phi_j(1, :))*(1/fs);
            phi_j(i, :) = (g_i / sqrt(Eg));
        end
    end

    % find remaining coefficients
    for i=1:M
        for j=i+1:M
            s_ij(i,j) = sum(s(i, :) .* phi_j(j, :)) * (1/fs);
        end
    end

    % return phis and coefficients
    x1(:, :) = phi_j; 
    x2(:, :) = s_ij;
end