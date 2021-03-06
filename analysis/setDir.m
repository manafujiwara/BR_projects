%% setDir
%%% First, cd to ~/scrip/analysis

addpath(genpath('./aux_files')) % chronax, eeglab are under the folder, Add the path when you need.
addpath(genpath('./decode/'))

% DIR.Toolbox = '/Applications/Psychtoolbox';
% addpath(genpath(DIR.Toolbox))

DIR.dataRaw       = '../../dataRaw';
DIR.cfg           = '../../config';
DIR.collectedData = '../../matRaw';
DIR.interpolated  = '../../interpolated';
DIR.epoched       = '../../epoched';

DIR.figures = '../../figures';
DIR.allSub  = '../../allSub';
DIR.Rdata   = '../../R/BR/results_eachTp';

DIR.figSlowPhase      = [DIR.figures '/slowPhase'];
DIR.figIMSortBPSlwPhs = [DIR.figures '/IM_sort_BP_a_SlwPhs'];
DIR.figIMBP           = [DIR.figures '/IM_BP'];
DIR.figRaw            = [DIR.figures '/IM_raw/'];
DIR.figGrpCmp         = [DIR.figures '/IM_grpCmp/'];
DIR.figRevision       = [DIR.figures '/IM_revision/'];


% NOT NEEDED ANYMORE --------------------------------
% DIR.filtGaze      = '../../filtredGaze';
% DIR.SVM     = '../../svm';
% DIR.velocity      = '../../velocity';
% DIR.slowPhase     = '../../slowPhase';
% DIR.figIMAUC          = [DIR.figures '/IM_AUC/'];
% DIR.figBRPerceptDrag  = [DIR.figures '/IM_BRPerceptDrag/'];
% DIR.figVilidity       = [DIR.figures '/validity'];
% DIR.figVelocity       = [DIR.figures '/velocity'  ];
% DIR.figIMSlwPhs       = [DIR.figures '/IM_SlwPhs'];
% DIR.figCentreGaze     = [DIR.figures '/IM_centreGaze/'];
% ---------------------------------------------------

% For continuous condition
DIR.excel = '../../epoched_excelForm';
DIR.clnData   = '../../CT_cleanData';
DIR.figCTGrpCmp = [DIR.figures '/CT_grpCmp/'];

pathCell = regexp(path, pathsep, 'split');
directory = fieldnames(DIR);

for idir = 1:length(directory)
    dirname = getfield(DIR, directory{idir});
    
    if ispc  % Windows is not case-sensitive
        onPath = any(strcmpi(dirname, pathCell));
    else
        onPath = exist(dirname, 'dir');
    end
    
    if onPath
        continue
    else
        if ~exist(dirname, 'dir')
            mkdir(dirname)
        end
        addpath(dirname)
    end
end


% DIR.clnData       = '../../CN_clnData/';
% DIR.figDecode         = [DIR.figures '/CN_decode'];
% DIR.figDecodeQual     = [DIR.figures '/CN_decodeQuality'];
% DIR.figAppearance     = [DIR.figures '/CN_appearance/'];

