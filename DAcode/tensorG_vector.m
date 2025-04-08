function G = tensorG_vector(k0,r1x,r1y,r1z,r2x,r2y,r2z)
% Green function for two-dimensional periodic array
eps = 1e30;
    dx = r2x - r1x;
    dy = r2y - r1y;
    dz = r2z - r1z;
   R12 = sqrt(dx.^2 + dy.^2 + dz.^2);
R12(R12 == 0) = eps;  % 防止除以零
    erx = dx ./ R12;
    ery = dy ./ R12;
    erz = dz ./ R12;

    common_term1 = (1 ./ R12 + 1i ./ (k0 * R12.^2) - 1 ./ (k0^2 * R12.^3));
    common_term2 = (-1 ./ R12 - 1i * 3 ./ (k0 * R12.^2) + 3 ./ (k0^2 * R12.^3));
    exp_term = exp(1i * k0 * R12) / (4 * pi);

    % 初始化 G 矩阵 (3 x 3 x N x N)
    G = zeros(3, 3, size(R12, 1), size(R12, 2));
    
    % 对角线部分
    % for idx = 1:3
    %     G(idx, idx, :, :) = common_term1 .* exp_term;
    % end
diag_mask = eye(3); % 对角掩码
G = diag_mask .* reshape(common_term1 .* exp_term, 1, 1, size(R12, 1), size(R12, 2));

    % 非对角线部分 (二阶项)
    G(1, 1, :, :) = squeeze(G(1, 1, :, :)) + common_term2 .* (erx .* erx) .* exp_term;
    G(1, 2, :, :) = common_term2 .* (erx .* ery) .* exp_term;
    G(1, 3, :, :) = common_term2 .* (erx .* erz) .* exp_term;

    G(2, 1, :, :) = common_term2 .* (ery .* erx) .* exp_term;
    G(2, 2, :, :) = squeeze(G(2, 2, :, :)) + common_term2 .* (ery .* ery) .* exp_term;
    G(2, 3, :, :) = common_term2 .* (ery .* erz) .* exp_term;

    G(3, 1, :, :) = common_term2 .* (erz .* erx) .* exp_term;
    G(3, 2, :, :) = common_term2 .* (erz .* ery) .* exp_term;
    G(3, 3, :, :) = squeeze(G(3, 3, :, :)) + common_term2 .* (erz .* erz) .* exp_term;

end


