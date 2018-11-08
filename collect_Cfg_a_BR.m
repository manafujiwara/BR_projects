function collect_Cfg_a_BR(DIR, subID)

for isub = 1 : length(subID)
    DIR.sub = [DIR.dataRaw '/' subID{isub} ];
    % Loading each subject once at a time,if you don't specify any return
    % variable to 'load' you'll just get 'Exp' in the workspace
    matFiles = dir(  [DIR.sub '/' subID{isub} '_*.mat' ] );
    matFiles = matFiles(find( ~cellfun(@isempty, regexp( {matFiles(:).name} , [subID{isub} '_(\d*)(\.mat)' ] , 'once' ) ) ) );
    
    imatFiles = 1;
    
    load([DIR.sub '/' matFiles(imatFiles).name])
    Cfg = Exp.Cfg;
    Cfg.samplingRate = 300;
    
    BR = Exp.BR;
    
    save([DIR.cfg '/'  subID{isub} '_Cfg_a_BR.mat'], 'Cfg', 'BR')
end