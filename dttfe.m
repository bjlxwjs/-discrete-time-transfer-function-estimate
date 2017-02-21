%% Copyright (C) 2009-2015   Lukas F. Reichlin
%%
%% This file is part of LTI Syncope.
%%
%% LTI Syncope is free software: you can redistribute it and/or modify
%% it under the terms of the GNU General Public License as published by
%% the Free Software Foundation, either version 3 of the License, or
%% (at your option) any later version.
%%
%% LTI Syncope is distributed in the hope that it will be useful,
%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%% GNU General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with LTI Syncope.  If not, see <http://www.gnu.org/licenses/>.
function [sys] = dttfe (dat, varargin)

  %% TODO: delays
  %% p: outputs,  m: inputs,  ex: experiments
  [~, p, m, ex] = size (dat);           % dataset dimensions

  varargin = horzcat (varargin(2:end), {'na'}, varargin(1), {'nb'}, varargin(1));

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
  Y = dat.y;
  Y=mat2cell(Y,(length(Y)));
  U = dat.u;
  U=mat2cell(U,(length(U)));
  tsam = dat.Ts;

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
  end

    A=A';
    B=B';
  sys = tf (B, A, tsam);                              % filt creates a transfer function in z^-1
end
