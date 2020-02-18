function seconds=converttoseconds(dateinput)

[~, ~, ~, H, MN, S] = datevec(dateinput);
seconds=H*3600+MN*60+S;

end