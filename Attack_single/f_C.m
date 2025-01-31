function [C,Vn] = f_C(n, measure_ind)
    C = zeros(n,n);
    for i = 1:size(measure_ind,2)
        C(measure_ind(i),measure_ind(i)) = 1;
    end
    Vn = C+eye(n)*0.00001;
  
end