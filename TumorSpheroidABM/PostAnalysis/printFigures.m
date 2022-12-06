function printFigures(f,cohort_name,varargin)

if nargin>2 && isstring(varargin{1})
    dir_path = varargin{1};
else
    dir_path = sprintf("../figs/%s",cohort_name);
end

if nargin>3
    reprint = varargin{2};
else
    reprint = false;
end

directories_made = false; % only attempt to make the directories once
file_formats = ["fig","png"];
for i = 1:numel(f)
    if ishandle(f(i)) && ~isempty(f(i).Name)
        if ~directories_made && ~exist("../figs","dir")
            mkdir("../figs")
        end
        for ffi = 1:length(file_formats)
            if ~directories_made && ~exist(sprintf("%s/%s",dir_path,file_formats(ffi)),"dir")
                mkdir(sprintf("%s/%s",dir_path,file_formats(ffi)))
            end
            if reprint || ~exist(sprintf("%s/%s/%s.%s",dir_path,file_formats(ffi),f(i).Name,file_formats(ffi)),"file")
                if file_formats(ffi)=="fig"
                    savefig(f(i),sprintf("%s/fig/%s",dir_path,f(i).Name))
                else
                    print(f(i),sprintf("%s/%s/%s",dir_path,file_formats(ffi),f(i).Name),sprintf("-d%s",file_formats(ffi)))
                end
            end
        end
        directories_made = true;
    end
end