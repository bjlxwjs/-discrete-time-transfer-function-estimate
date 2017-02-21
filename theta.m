function Theta = theta (Phi, Y, i, n)
    
  if (numel (Phi) == 1)                             % single-experiment dataset
    %% use 'square-root algorithm'
    A = horzcat (Phi{1}, Y{1}(n(i)+1:end, i));      % [Phi, Y]
    R0 = triu (qr (A, 0));                          % 0 for economy-size R (without zero rows)
    R1 = R0(1:end-1, 1:end-1);                      % R1 is triangular - can we exploit this in R1\R2?
    R2 = R0(1:end-1, end);
    Theta = ls_svd (R1, R2);                    % R1 \ R2
    
    %% Theta = Phi \ Y(n+1:end, :);                 % naive formula
    %% Theta = __ls_svd__ (Phi{1}, Y{1}(n(i)+1:end, i));
  else                                              % multi-experiment dataset
    %% TODO: find more sophisticated formula than
    %% Theta = (Phi1' Phi1 + Phi2' Phi2 + ...) \ (Phi1' Y1 + Phi2' Y2 + ...)
    
    %% covariance matrix C = (Phi1' Phi + Phi2' Phi2 + ...)
    tmp = cellfun (@(Phi) Phi.' * Phi, Phi, 'uniformoutput', false);
    %% rc = cellfun (@rcond, tmp);                     % also test  QR or SVD
%     C = plus (tmp{:});
      C =  (tmp{1});

    %% PhiTY = (Phi1' Y1 + Phi2' Y2 + ...)
%     tmp = cellfun (@(Phi, Y) Phi.' * Y(n(i)+1:end, i), Phi, Y, 'uniformoutput', false);
    tmp = cellfun (@(Phi) Phi.' * Y(n(i)+1:end, i), Phi, Y, 'uniformoutput', false);
    PhiTY = (tmp{1});
    
    %% pseudoinverse  Theta = C \ Phi'Y
    Theta = ls_svd(C, PhiTY);
  end
  
end
