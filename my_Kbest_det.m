format short g
% init the parameter
H = [
(376    +26j) , ( 32     -88j);
(16    -34j) , (  444 -   20j);
];

order = 4;

Data_in = [
    (155   -349j);
    (-365 +   186j);
];

K_best = 8;

% start

y_zf = pinv(H) * Data_in;

switch order
  case 2 % QPSK
    radius = norm((Data_in - y_zf) .* sqrt(2));
    Data_in = Data_in .* sqrt(2);
    xMax = 1 + 1j;
  case 4 % 16QAM
    radius = norm((Data_in - y_zf) .* sqrt(10)); 
    Data_in = Data_in .* sqrt(10);
    xMax = 3 + 3j;
  case 6 % 64QAM
    radius = norm((Data_in - y_zf) .* sqrt(42)); 
    Data_in = Data_in .* sqrt(42);
    xMax = 7 + 7j;
end

[Out, wHat] = Kbest_det(2, xMax, Data_in, H, radius, K_best);

function [xHat, wHat] = Kbest_det(m, xMax, Data_in, H, radius, k_value)

  % n: antenna
  % m: layer
  % Complex -> real 

  m = 2 * m;
  H = [[real(H), -imag(H)]; 
       [imag(H),  real(H)]];
  Data_in = [real(Data_in); imag(Data_in)];  
  xMax = [real(xMax); imag(xMax)];
  xHat = -xMax(1) : 2 : xMax(1); % range in QAM candidate for each QAM

  [Hq, Hr] = qr(H);
  yR = Hq' * Data_in;

  dp = radius.^2;
  level = 1;
  dim = m;
  dimp = dim + 1;

  
  dp  = dp .* ones(1, size(xHat, 2)); % radius 1x4
  
  py = yR;
  py_extend = py(:, ones(1, size(xHat, 2))); % 4x1 to 4x4
  py = py_extend;

  for numLayer = 1 : m
    k = k_value; 
    if (length(dp) < k)
      k = length(dp);
    end
    
    xHat_next = [];
    py_next = [];
    dp_next = [];

    dim = dim - 1;
    dimp = dimp - 1;
    
    tmp_w = (py(dimp,:) .* ones(1, size(xHat,2)) - Hr(dimp, dimp) .* xHat(level, :)).^2;
    dp = dp - tmp_w;

    for i = 1 : length(dp)
      if (dp(i) < 0)
        dp(i) = 0;
      end
    end

    if (length(dp) < k)
      sort_idx = 1:length(dp);
    else
      [~, sort_idx] = sort(dp);
      sort_idx = sort_idx(length(dp) - k + 1: length(dp));
    end

    xHat = xHat(:, sort_idx);
    py = py(:, sort_idx);
    dp = dp(sort_idx);

    if (numLayer == m)
      if (size(xHat, 2) == 0)
        xHat = Data_in \ H; 
        wHat = 0;
      else
        inverse_xHat = xHat(:, find(dp == max(dp)));
        wHat = min(dp);
        xHat = zeros(m,1);
        for ii = 1:m
          xHat(ii) = inverse_xHat(m+1-ii);
        end
      end
      break;
    end

    % for the next iteration
    for ii = 1 : size(xHat, 2)
      tmp_py = py(1:dim, ii) - Hr(1:dim,dimp) * xHat(level, ii);  % the residue target

      tmp_range = -xMax(1):2:xMax(1);
      
      tmp_xHat = xHat(1:level, ii) * ones(1, length(tmp_range));
      tmp_xHat = [tmp_xHat; tmp_range];
      xHat_next = [xHat_next tmp_xHat];

      dp_next  = [dp_next ones(1, length(tmp_range)) .* dp(ii)];

      py_extend = tmp_py(1:(m - level), 1) * ones(1, size(tmp_range, 2));
      py_next = [py_next py_extend];

    end

    dp = dp_next;  % raduis
    py = py_next;  % yR
    xHat = xHat_next; % record for final answer
    level = level + 1;
  end
  m = m / 2;
  xHat = xHat(1:m)+1i*xHat(m+1:2*m);
end
