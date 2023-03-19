function aclr = manualSetIgnore(aclr)


%MANUALSETIGNORE Manual update ignore flag
%   Prasanna Sritharan, August 2020


% user settings
addpath('..');
user = getUserScriptSettings();
outpath = user.OUTPATH;


% set ignore (0=include,1=ignore)
aclr.ACLR.SUBJECT17.HOPP02.ignore = 1;
aclr.ACLR.SUBJECT17.HOPP04.ignore = 1;
aclr.ACLR.SUBJECT21.HOPP02.ignore = 1;
aclr.ACLR.SUBJECT21.HOPP04.ignore = 1;
aclr.ACLR.SUBJECT26.HOPP05.ignore = 1;

% save metastruct
save([outpath '\aclr.mat'],'aclr');

end

