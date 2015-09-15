function Matlab2CSV4R(header,SavePrompt,varargin)
% converts data into a CSV file suitable for input into R
%   varargin are numeric vectors or string arrays which must have the same length as header
%   header is a string array with the same length as the number of input
%   vectors
n=nargin;
if length(header)~=(n-2)
    error('Number of vectors inconsistent with header');
end
l=length(varargin{1});
for i=2:n-2;
    if length(varargin{i})~=l
        error('Length of input vectors not the same');
    end
end
%get output file
[FileName,PathName,FilterIndex] = uiputfile('*.txt',SavePrompt);
fid=fopen([PathName FileName],'w');
% create text file
%write header
fprintf(fid,cat(2,repmat('\t%s',[1,n-2]),'\n\r'),header{:});
%determine format string
fs='%u'; % row number
for i=1:n-2
    if iscell(varargin{i})
        fs=cat(2,fs,'\t%s');
    else
        fs=cat(2,fs,'\t%10.8e');
    end
end
fs=cat(2,fs,'\n\r');
for i=1:l
    linedata{1}=i;
    for j=1:n-2
        if isa(varargin{j},'cell')
            linedata{1+j}=varargin{j}{i};
        else
            linedata{1+j}=varargin{j}(i);
        end
    end
    fprintf(fid,fs,linedata{:});
    if mod(i,100)==0
        display(['Line ' num2str(i) ' of ' num2str(l)]);
    end
end      
fclose(fid);
end

