function [W, pi,sigma] = ReadDataFile(path)
 if ~exist('path','var')
   file1 = '/Users/shuchenwu/Downloads/bsc-et-2/examples/W128.csv';
 else 
   file1 = path;
 end



fid = fopen(file1);

tline = fgets(fid);
sigma = textscan(tline,'%f');
sigma = sigma{1};
tline = fgets(fid);

pi = textscan(tline,'%f');
pi = pi{1};

W = [];
tline = fgets(fid);

while (ischar(tline))

  scan = textscan(tline,'%f','Delimiter',',');
  W = [W;scan{1}'];
  tline = fgets(fid);
end
W = W';
fclose(fid);