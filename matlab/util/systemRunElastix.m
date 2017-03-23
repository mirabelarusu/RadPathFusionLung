function systemRunElastix(elastix_exe, fIm, mIm, t0, param, path_out)
    % create folder for output if it doesn't exists
    if (exist(path_out)==0)
        disp(['Creating folder', path_out]);
        mkdir(path_out)
    else
        disp(['Writing to folder: "', path_out,'"']);
    end
    

    el_cmd = [elastix_exe,' -p ', param,' -f ',fIm, ' -m ', mIm, ' -out ', path_out];
    if isempty(t0)==false
        el_cmd = [el_cmd,' -t0 ',t0];
    end
    
    el_cmd = [el_cmd,' >/dev/null'];
    system(el_cmd);

end