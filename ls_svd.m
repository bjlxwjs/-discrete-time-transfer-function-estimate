function x = ls_svd (A, b)

  %% solve the problem Ax=b
  %% x = A\b  would also work,
  %% but this way we have better control and warnings
  %% solve linear least squares problem by pseudoinverse
  %% the pseudoinverse is computed by singular value decomposition
  %% M = U S V*  --->  M+ = V S+ U*
  %% Th = Ph \ Y = Ph+ Y
  %% Th = V S+ U* Y,   S+ = 1 ./ diag (S)

  [U, S, V] = svd (A, 0);                           % 0 for 'economy size' decomposition
  S = diag (S);                                     % extract main diagonal
  r = sum (S > eps*S(1));
  if (r < length (S))
    warning ('arx: rank-deficient coefficient matrix');
    warning ('sampling time too small');
    warning ('persistence of excitation');
  end
  V = V(:, 1:r);
  S = S(1:r);
  U = U(:, 1:r);
  x = V * (S .\ (U' * b));                          % U' is the conjugate transpose

  end