function [sys] = dttfe (dat, varargin)

  %% TODO: delays
  %% p: outputs,  m: inputs,  ex: experiments
  [~, p, m, ex] = size (dat);           % dataset dimensions

%   if (is_real_scalar (varargin{1}))     % arx (dat, n, ...)
     varargin = horzcat (varargin(2:end), {'na'}, varargin(1), {'nb'}, varargin(1));
%   end

  if (isstruct (varargin{1}))           % arx (dat, opt, ...), arx (dat, n, opt, ...)
    varargin = horzcat (opt2cell (varargin{1}), varargin(2:end));
  end

  nkv = numel (varargin);               % number of keys and values
  

  %% default arguments
  na = [];
  nb = [];
  nk = 0;

  %% handle keys and values
  for k = 1 : 2 : nkv
    key = lower (varargin{k});
    val = varargin{k+1};
    switch (key)
      case 'na'
        na = val;
      case 'nb'
        nb = val;
    end
  end

  %% extract data  
  %%读取数据并转换为cell格式
  Y = dat.y;
  Y=mat2cell(Y,(length(Y)));
  U = dat.u;
  U=mat2cell(U,(length(U)));
  tsam = dat.Ts;
%   tsam=mat2cell(tsam,(length(tsam)));

  %% multi-experiment data requires equal sampling times  
  
   % tsam = tsam{1};

  
    na = repmat (na, p, 1);                         % na(p-by-1)
    nb = repmat (nb, p, m);                         % nb(p-by-m)
  

  max_nb = max (nb, [], 2);                         % one maximum for each row/output, max_nb(p-by-1)
  n = max (na, max_nb);                             % n(p-by-1)

  %% create empty cells for numerator and denominator polynomials
  num = cell (p, m+p);
  den = cell (p, m+p);

  %% MIMO (p-by-m) models are identified as p MISO (1-by-m) models
  %% For multi-experiment data, minimize the trace of the error
 i=1;                              % for every output
    Phi = cell (ex, 1);                             % one regression matrix per experiment
    for e = 1 : ex                                  % for every experiment  
      %% avoid warning: toeplitz: column wins anti-diagonal conflict
      %% therefore set first row element equal to y(1)
      PhiY = toeplitz (Y{e}(1:end-1, i), [Y{e}(1, i); zeros(na(i)-1, 1)]);
      %% create MISO Phi for every experiment
      PhiU = arrayfun (@(x) toeplitz (U{e}(1:end-1, x), [U{e}(1, x); zeros(nb(i,x)-1, 1)]), 1:m, 'uniformoutput', false);
      PhiU=cell2mat(PhiU);
      Phi = (horzcat (-PhiY, PhiU));
      Phi(1:n(i)-1,:)=[];
      Phi=mat2cell(Phi,(length(Phi)));
    end
 for i = 1 : p       
    %% compute parameter vector Theta
    Theta = theta (Phi, Y, i, n);

    %% extract polynomial matrices A and B from Theta
    %% A is a scalar polynomial for output i, i=1:p
    %% B is polynomial row vector (1-by-m) for output i
    A = [1; Theta(1:na(i))];                                % a0 = 1, a1 = Theta(1), an = Theta(n)
    ThetaB = Theta(na(i)+1:end);                            % all polynomials from B are in one column vector
    B=ThetaB;

%     B = mat2cell (ThetaB, nb(i,:));                         % now separate the polynomials, one for each input
%     B = reshape (B, 1, []);                                 % make B a row cell (1-by-m)
%     B = cellfun (@(B) [zeros(1+nk, 1); B], B, 'uniformoutput', false);  % b0 = 0 (leading zero required by filt)
%     B=cell2mat(B);
    
    %% add error inputs
%     Be = repmat ({0}, 1, p);                                % there are as many error inputs as system outputs (p)
%     Be(i) = [zeros(1,nk), 1];                               % inputs m+1:m+p are zero, except m+i which is one
%     num(i, :) = [B, Be];                                    % numerator polynomials for output i, individual for each input
%     den(i, :) = repmat ({A}, 1, m+p);                       % in a row (output i), all inputs have the same denominator polynomial
  end

  %% A(q) y(t) = B(q) u(t) + e(t)
  %% there is only one A per row
  %% B(z) and A(z) are a Matrix Fraction Description (MFD)
  %% y = A^-1(q) B(q) u(t) + A^-1(q) e(t)
  %% since A(q) is a diagonal polynomial matrix, its inverse is trivial:
  %% the corresponding transfer function has common row denominators.
    A=A';
    B=B';
  sys = tf (B, A, tsam);                              % filt creates a transfer function in z^-1

  %% compute initial state vector x0 if requested
  %% this makes only sense for state-space models, therefore convert TF to SS
  
%   if (nargout > 1)
%     sys = prescale (ss (sys(:,1:m)));
%     x0 = __sl_ib01cd__ (Y, U, sys.a, sys.b, sys.c, sys.d, 0.0);
%     %% return x0 as vector for single-experiment data
%     %% instead of a cell containing one vector
%     if (numel (x0) == 1)
%       x0 = x0{1};
%     endif
%     varargout{1} = x0;
%   endif

end
