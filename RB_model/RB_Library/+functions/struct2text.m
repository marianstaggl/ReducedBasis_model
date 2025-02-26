function struct2text(stru, filename)
%STRUCT2TEXT write structdata into a textfile
% write in try and catch block
try
    fid = fopen(filename, 'w');
    write_struct(fid, stru); fclose(fid);
catch me
    fclose(fid);
    rethrow(me)
end
end

function write_struct(fid, stru)
% get the fields of the struct
s_f = fields(stru);
for i=1:numel(s_f)
    stri = [s_f{i}, ': ' num2str(stru.(s_f{i})) '\n'];
    fprintf(fid, stri);
end
end