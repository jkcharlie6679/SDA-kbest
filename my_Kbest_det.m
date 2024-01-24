format short g
H = [
(334    -30j) , ( 40     -6j);
( -6     -4j) , (  378 +    0j);
];

order = 4;

Data_in = [
    (591   -544j);
    (-223 +   93j);
];

if (order == 6)
  Data_in = Data_in ./ 12;
end

y_zf = adjoint(H) * Data_in / 256;

switch order
  case 2
    radius = norm((Data_in - y_zf).*362/256); % sqrt(2) nrom -> sqrt(|a1| + |a2|)
    Data_in = Data_in.* 364 / 256;% sqrt(2)
    xMax = 1 + 1j;
  case 4
    radius = norm((Data_in - y_zf).*810) / 256; % sqrt(10)
    Data_in = Data_in.* 810 / 256; % sqrt(10)
    xMax = 3 + 3j;
  case 6
    radius = norm((Data_in - y_zf).* 1659) / 256; % sqrt(42)
    Data_in = Data_in.* 1659 / 256; % sqrt(42)
    xMax = 7 + 7j;
end

[Out, wHat] = Kbest_det(2, xMax, Data_in, H, radius, 4);

function [xHat, wHat] = Kbest_det(m, xMax, Data_in, H, radius, k_value)

  % n: antenna
  % m: layer
  % Complex -> real 

  m = 2 * m;
  H = [[real(H), -imag(H)]; 
       [imag(H),  real(H)]];
  Data_in = [real(Data_in); imag(Data_in)];  
  xMax = [real(xMax); imag(xMax)];
  xHat = -xMax(1):2:xMax(1); % range in QAM candidate for each QAM

  [Hq, Hr] = qr(H);

  yR = Hq' * Data_in;

  % dp = radius.^2;
  level = 1;
  dim = m;
  dimp = dim + 1;

  py = yR;
  ped = zeros(1, length(xHat));
  % dp  = dp .* ones(1, size(xHat, 2)); % radius 1x4

  py_extend = py(:, ones(1, size(xHat, 2))); % 4x1 to 4x4
  py = py_extend;


  for numLayer = 1 : m
    k = k_value; 
    if (length(ped) < k)
      k = length(ped);
    end
    
    xHat_next = [];
    py_next = [];
    ped_next = [];
    % dp_next = [];

    dim = dim - 1;
    dimp = dimp - 1;
    
    tmp_w = (py(dimp,:) .* ones(1, size(xHat,2)) - Hr(dimp, dimp) .* xHat(level, :)).^2;
    % tmp_dp = (dp) - tmp_w; % for next radius 
    
    ped = ped + tmp_w;

    if (length(ped) < k)
      sort_idx = 1:length(ped);
    else
      [~, sort_idx] = sort(ped);
      sort_idx = sort_idx(1:k);
    end

    xHat = xHat(:, sort_idx);
    py = py(:, sort_idx);
    % tmp_dp = tmp_dp(sort_idx);
    ped = ped(sort_idx);
    if (numLayer == m)
      if (size(xHat, 2) == 0)
        xHat = Data_in\H; 
        wHat = 0;
      else
        inverse_xHat = xHat(:, find(ped==max(ped)));
        wHat = min(ped);
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

      ped_next = [ped_next ones(1,length(tmp_range)) .* ped(ii)];

      % dp_next  = [dp_next ones(1, length(tmp_range)) .* tmp_dp(ii)];

      py_extend = tmp_py(1:(m - level), 1) * ones(1, size(tmp_range, 2));
      py_next = [py_next py_extend];

    end

    % dp = dp_next;  % raduis
    py = py_next;  % yR
    ped = ped_next; % add for sort
    xHat = xHat_next; % record for final answer

    level = level + 1;
  end
  m = m/2;
  xHat = xHat(1:m)+1i*xHat(m+1:2*m);
end
