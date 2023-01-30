Load data
Train RF with OOB estimators
Train RF without OOB estimators
Fitting models
Save AMP model
Predict test
# AMP (unbalanced training) CLASSIFIER 
Accuracy: 0.926
MCC: 0.862

Confusion matrix:
[[784 136]
 [  0 920]]

              precision    recall  f1-score   support

         AMP      1.000     0.852     0.920       920
        NAMP      0.871     1.000     0.931       920

    accuracy                          0.926      1840
   macro avg      0.936     0.926     0.926      1840
weighted avg      0.936     0.926     0.926      1840



Load bench set
Train bench set
Test bench set
# AMP (benchmark training) CLASSIFIER 
Accuracy: 0.946
MCC: 0.893

Confusion matrix:
[[846  74]
 [ 26 894]]

              precision    recall  f1-score   support

         AMP      0.970     0.920     0.944       920
        NAMP      0.924     0.972     0.947       920

    accuracy                          0.946      1840
   macro avg      0.947     0.946     0.946      1840
weighted avg      0.947     0.946     0.946      1840



Load datasets
Fit Hemo sets
Save Hemo model
# Hemo CLASSIFIER 
Accuracy: 0.936
MCC: 0.873

Confusion matrix:
[[102   8]
 [  6 104]]

              precision    recall  f1-score   support

        Hemo      0.944     0.927     0.936       110
     NonHemo      0.929     0.945     0.937       110

    accuracy                          0.936       220
   macro avg      0.937     0.936     0.936       220
weighted avg      0.937     0.936     0.936       220

January 30th 2023
scikit-learn 1.1.2
