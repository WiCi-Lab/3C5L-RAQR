function key = make_eit_cache_key(tag, vals)
vals = vals(:).';
key = [char(tag), '|', sprintf('%.17g,', vals)];
end