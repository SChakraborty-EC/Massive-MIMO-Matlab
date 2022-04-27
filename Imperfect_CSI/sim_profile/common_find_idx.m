function idxhat = common_find_idx(ML_symbol, symbol_list)

idxhat = zeros(length(ML_symbol),1);

for i = 1 : length(ML_symbol)
    idxhat(i,1) =  find(symbol_list == ML_symbol(i,1));
end
