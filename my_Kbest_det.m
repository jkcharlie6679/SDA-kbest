format short g
H = [
( 360 - 0j), ( 4 +   0j);
(  2 - 4j), ( 358 -  4j);
];

order = 6;

Data_in = [
    (572  -3123j);
    (-3145 + 4429j);
];

if (order == 6)
  Data_in = Data_in ./ 12;
end
% Data input from OAI finish

% ZF should be pinv(H) * Data_in

y_zf = pinv(H) * Data_in * det(H) / 256;

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
% disp(radius)

% disp(Data_in);
[Out, wHat] = Kbest_det(2, xMax, Data_in, H, radius, 4);

function [xHat, wHat] = Kbest_det(m, xMax, Data_in, H, radius, k_value)

  % n: antenna
  % m: layer

  % xMax = xMax * ones(m,1);

  % Complex -> real conversion if necessary

  m = 2 * m;
  H = [[real(H), -imag(H)]; 
       [imag(H),  real(H)]];
  Data_in = [real(Data_in); imag(Data_in)];  
  xMax = [real(xMax); imag(xMax)];
  xHat = -xMax(1):2:xMax(1); % range in QAM candidate for each QAM

  [Hq,Hr] = qr(H);

  yR      = Hq'*Data_in;

  % disp(Hq)
  % disp(yR)
  dp = radius;
  level = 1;
  dim = m;
  dimp = dim + 1;

  py = yR;

  

  ped = zeros(1, length(xHat));
  dp  = dp .* ones(1, size(xHat, 2)); % radius 1x4

  py_extend = py(:, ones(1, size(xHat, 2))); % 4x1 to 4x4
  py = py_extend;


  for numLayer = 1:m
    k = k_value; 
    if (length(ped) < k)
      k = length(ped);
    end
    
    xHat_next = [];
    py_next = [];
    ped_next = [];
    dp_next = [];

    dim = dim - 1;
    dimp = dimp - 1;
    
    tmp_w = (py(dimp,:) .* ones(1, size(xHat,2)) - Hr(dimp, dimp) .* xHat(level, :)) .^ 2 ;
    % disp(tmp_w)
    tmp_dp = (dp).^2 - tmp_w; % 差距
    
    ped = ped + tmp_w; % P_k
    % disp(ped)
    if (length(ped) < k)
      sort_idx = 1:length(ped);
    else
      [sort_ped, sort_idx] = sort(ped);
      disp(sort_idx)
      sort_idx = sort_idx(1:k);
      % sort_idx = sort(sort_idx);
    end
    % if (numLayer == 1)
        
    % end
    tmp_dp = sqrt(tmp_dp);
    
    % if (numLayer > 0)
    %     fprintf("This is layer %d\n", numLayer)
    %     disp("xHat: ")
    %     disp(xHat);
    %     disp("py: ")
    %     disp(py);
    %     disp("tmp_dp: ")
    %     disp(tmp_dp);
    %     disp("ped: ")
    %     disp(ped);
    % end
    disp(xHat)
    xHat = xHat(:,sort_idx);
    py = py(:,sort_idx);
    tmp_w = tmp_w(sort_idx);
    tmp_dp = tmp_dp(sort_idx);
    ped = ped(sort_idx);
    % if (numLayer > 0)
    %     disp("Afer sorting: ");
    %     disp("xHat: ")
    %     disp(xHat);
    %     disp("py: ")
    %     disp(py);
    %     disp("tmp_dp: ")
    %     disp(tmp_dp);
    %     disp("ped: ")
    %     disp(ped);
    % end

    if (numLayer == m)
        
      if (size(xHat,2)==0)
        xHat = inv(H) * Data_in; 
        wHat = 0;
      else
        inverse_xHat = xHat(:,find(ped==min(ped)));
        wHat = min(ped);
        xHat = zeros(m,1);
        for ii = 1:m
          xHat(ii) = inverse_xHat(m+1-ii);
        end
        
      end
      break;
    end

    % disp("--------------------------------------");
    for ii = 1:size(xHat, 2)
      % fprintf("---This is ii: %d---\n", ii);
      tmp_py = py(1:dim,ii) - Hr(1:dim,dimp)*xHat(level,ii);  % the residue target
      % disp("tmp_py: ")
      % disp(tmp_py)
      tmp_range = -xMax(1):2:xMax(1);
      
      tmp_xHat = xHat(1:level,ii)*ones(1,length(tmp_range));
      % disp("tmp_xHat: ")
      % disp(tmp_xHat)
      tmp_xHat = [ tmp_xHat; tmp_range];
      % disp(tmp_xHat)
      xHat_next = [xHat_next tmp_xHat];
      % disp("xHat_next: ")
      % disp(xHat_next)

      ped_next = [ped_next ones(1,length(tmp_range)).*ped(ii)];
      % disp("ped_next: ")
      % disp(ped_next)

      dp_next  = [dp_next ones(1,length(tmp_range)).*tmp_dp(ii)];
      % disp("dp_next: ")
      % disp(dp_next)
      
      py_extend = tmp_py(1:(m-level),1)*ones(1,size(tmp_range,2));
      py_next = [py_next py_extend];
      % disp("py_extend: ")
      % disp(py_next)
    end

    dp = dp_next;
    py = py_next;
    ped = ped_next;
    xHat = xHat_next;

    level = level + 1;
    disp("--------------------------------------");
    disp(ped)
  end
  

  m = m/2;
  xHat = xHat(1:m)+1i*xHat(m+1:2*m);
end
