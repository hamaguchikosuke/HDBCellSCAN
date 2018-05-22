function [output_dir,SubDirs] = get_output_dir_from_db(db,sep)
% [output_dir,SubDirs] = get_output_dir_from_db(db,sep)
%  from db (database file for imaging analysis use for Suite2P)
%  generate output directory name.
%  sep: separater of output_dir name.
% ex) 
% 
%  KH 20170909

if  ischar(db.expts)
        expts=db.expts;
        db=rmfield(db,'expts')
        db.expts{1}   = expts;
end

SubDirs=[];
for k = 1:length(db.expts)
    if iscell(db.expts)
        SubDirs{k}   = db.expts{k};
    elseif isvector(db.expts)
        SubDirs{k}    = num2str(db.expts(k));
    end
end


CharSubDirs = '';
for i = 1:length(SubDirs)
    CharSubDirs = [CharSubDirs SubDirs{i} sep];
end
CharSubDirs = CharSubDirs(1:end-length(sep));

output_dir = CharSubDirs;
