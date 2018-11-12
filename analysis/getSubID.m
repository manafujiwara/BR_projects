folders_sub = dir( [ DIR.interpolated '/'  ] );
% folders_sub = dir( [ DIR.dataRaw '/'  ] );

expMatch = '<\[cp]0\d\d\d';
loadFileIdx = find( ~cellfun( @isempty, regexp( {folders_sub(:).name} , expMatch, 'match' ) ) );


%% get subID + remove zero from file name
if isempty(loadFileIdx) % zero has been removed from all the file names
    
    expMatch = '\<[cp]\d{3,3}';
    loadFileIdx = find( ~cellfun( @isempty, regexp( {folders_sub(:).name} , expMatch, 'match' ) ) );
    tmp = char( folders_sub(loadFileIdx).name ) ;
    subID = cellstr( tmp(:,1:4) );

else % there is/are a file/files just added with zero
    subID_ori = char( folders_sub(loadFileIdx).name );
    subID = cellstr([ subID_ori(:,[1,3:5]) ]);
    subID_ori = cellstr(subID_ori);
    
    for iFolder = 1: length(subID)
        
        DIR.sub = [ DIR.dataRaw '/' subID{iFolder} ];
        movefile( [ DIR.dataRaw '/' subID_ori{iFolder}], DIR.sub )
        
        files_sub = dir( DIR.sub );
        
        expMatch = ['\<' subID_ori{iFolder}];
        
        loadFileIdx = find( ~cellfun( @isempty, regexp( {files_sub(:).name} , expMatch, 'match' ) ) );
        files_ori = char(files_sub(loadFileIdx).name) ;
        files = cellstr( files_ori( :,[1,3:end] ) );
        files_ori = cellstr(files_ori);
        
        for iFile = 1: length(loadFileIdx)
            movefile( [ DIR.sub '/' files_ori{iFile} ], [ DIR.sub '/' files{iFile}] )
        end
        
    end
    
    expMatch = '\<[cp]\d{3,3}';
    loadFileIdx = find( ~cellfun( @isempty, regexp( {folders_sub(:).name} , expMatch, 'match' ) ) );
    subID_ori = char( folders_sub(loadFileIdx).name );
    subID = cellstr([ subID_ori(:,[1,3:5]) ]);
    
end