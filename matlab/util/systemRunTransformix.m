function systemRunTransformix(exe, mIm, outParam, path_out)
    % create folder for output if it doesn't exists
    if (exist(path_out)==0)
        disp(['Creating folder', path_out]);
        mkdir(path_out)
    else
        disp(['Writing to folder: "', path_out,'"']);
    end
    
    tr_cmd = [exe,' -tp ', outParam, ' -in ', mIm, ' -out ', path_out, ...
        ' >/dev/null'];
    system(tr_cmd);

end