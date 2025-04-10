function [F1_macro,F1_micro]=compute_f1_scores(Y, A)
%%% A (predicted labels) and Y (true labels)
% Convert categorical labels (1-15) into one-hot binary matrices
    num_classes = 15;
    num_subjects = numel(Y);
    
    % Convert A and Y to one-hot encoding
    Y_onehot = zeros(num_classes, num_subjects);
    A_onehot = zeros(num_classes, num_subjects);
    
    for j = 1:num_subjects
        Y_onehot(Y(j), j) = 1; % Mark true label in Y
        A_onehot(A(j), j) = 1; % Mark predicted label in A
    end

    % Ensure matrices are logical
    Y_onehot = logical(Y_onehot);
    A_onehot = logical(A_onehot);

    % Initialize macro-F1 components
    f1_macro_list = zeros(num_classes, 1);
    
    % Compute global TP, FP, FN for micro-F1
    TP_global = 0; FP_global = 0; FN_global = 0;

    % Compute class-wise F1 scores (for Macro-F1)
    for i = 1:num_classes
        TP = sum(Y_onehot(i, :) & A_onehot(i, :)); % True Positives
        FP = sum(~Y_onehot(i, :) & A_onehot(i, :)); % False Positives
        FN = sum(Y_onehot(i, :) & ~A_onehot(i, :)); % False Negatives

        % Compute precision and recall for class i
        P = TP / (TP + FP + eps); % Avoid division by zero
        R = TP / (TP + FN + eps);

        % Compute F1 for class i
        F1 = 2 * (P * R) / (P + R + eps);
        f1_macro_list(i) = F1;

        % Accumulate for micro-F1
        TP_global = TP_global + TP;
        FP_global = FP_global + FP;
        FN_global = FN_global + FN;
    end

    % Macro-F1 score (average across classes)
    F1_macro = mean(f1_macro_list);

    % Micro-F1 score (global precision and recall)
    P_micro = TP_global / (TP_global + FP_global + eps);
    R_micro = TP_global / (TP_global + FN_global + eps);
    F1_micro = 2 * (P_micro * R_micro) / (P_micro + R_micro + eps);
end
