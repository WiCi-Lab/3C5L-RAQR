function [model_cfg, args] = extract_model_cfg(args)
if ~isempty(args) && isstruct(args{end})
    model_cfg = args{end};
    args(end) = [];
else
    model_cfg = default_model_cfg();
end
end
