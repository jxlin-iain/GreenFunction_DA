function G = tensorG_mnp(u1, u2, N)
    % Computes a 3x3xNxN tensor G based on two input 3xN matrices u1 and u2
    %
    % Inputs:
    %   u1: 3xN matrix (each column represents a 3D vector for one point)
    %   u2: 3xN matrix (each column represents a 3D vector for one point)
    %
    % Output:
    %   G: 3x3xNxN tensor where G(:,:,i,j) = u1(:,i) * conj(u2(:,j))'

    % Ensure input dimensions are valid

    u1_expanded = reshape(u1,[3,1,N^2]);
    u2_expanded = reshape(conj(u2),[1,3,N^2]);
    G = bsxfun(@times,u1_expanded,u2_expanded);
end