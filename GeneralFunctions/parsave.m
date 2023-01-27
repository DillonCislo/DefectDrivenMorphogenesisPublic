function parsave(fname, varname, x)
% By Dillon Cislo

if iscell(varname)
    
    for i = 1:numel(varname)
        varname{i} = genvarname(varname{i});
        eval([varname{i} '= x{i};']);
    end
    varname = cellfun(@(x) ['"' x '", '], varname, 'Uni', false);
    varname = [varname{:}];
    varname((end-1):end) = [];
    eval(['save("' fname '", ' varname ');']);
        
else
    
    varname = genvarname(varname);
    eval([varname '= x;']);
    save(fname, varname);
    
end

end
