function f1_scores=Macro_F1_score(Ypred,Ytest)
% Macro-F1 score
num_classes=length(unique(Ytest));
f1s = zeros(num_classes, 1);

for c = 1:num_classes
    class = c;
    tp = sum(Ypred == class & Ytest == class);
    fp = sum((Ypred == class) & (Ytest ~= class));
    fn = sum((Ypred ~= class) & (Ytest == class));
    
    if tp + fp == 0 || tp + fn == 0
        f1 = 0;
    else
        precision = tp / (tp + fp);
        recall = tp / (tp + fn);
        if precision + recall == 0
            f1 = 0;
        else
            f1 = 2 * precision * recall / (precision + recall);
        end
    end
    f1s(c) = f1;
end
f1_scores = mean(f1s);