clear all;
clc;
fid = fopen('result.txt');
result = zeros(101,101);
no_of_rows = 101;
no_of_col = 101;
k = 1;
while k < 15000
    for i=k:(k + 100)
        for j=1:101
            result(i,j) = fscanf(fid,'%f',1);
        end
    end
    k = k*101 + 1;
end

	
	