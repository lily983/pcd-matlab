function v = skew2vec(S)
    % Extract vector from 3x3 skew-symmetric matrix
    % Input: 3x3 skew-symmetric matrix S
    % Output: 3x1 vector v
    
    % Ensure S is 3x3 skew-symmetric
    if ~isequal(size(S), [3 3]) || ~isequal(S, -S.')
        error('Input must be a 3x3 skew-symmetric matrix.');
    end

    % Extract vector components
    v = [S(3, 2); S(1, 3); S(2, 1)];
end
