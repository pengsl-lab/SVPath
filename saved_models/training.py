import numpy as np
import pandas as pd
import argparse
import prettytable as pt
from imblearn.over_sampling import BorderlineSMOTE
from sklearn.metrics import matthews_corrcoef, accuracy_score, confusion_matrix
from sklearn.model_selection import ShuffleSplit
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn import svm
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import GradientBoostingClassifier
import pickle
import os

# from sklearn import metrics
# import csv
# import vcf
# from sklearn.metrics import precision_score, recall_score, f1_score

def get_features(features_csv):
    data = pd.read_csv(features_csv)
    cols = data.columns.to_list()
    data = pd.DataFrame(data)
    data = data.values.tolist()
    num_pathogenic = 0
    num_benign = 0
    for i in range(len(data)):
        if data[i][-1] == 'Pathogenic':
            data[i][-1] = 1
            num_pathogenic += 1
        elif data[i][-1] == 'Benign':
            data[i][-1] = 0
            num_benign += 1
    locationAndId = []
    if 'Location' in cols and 'ID' in cols and 'SV_LEN' in cols:
        locationAndId = [x[:3] for x in data]
        X = [x[3:-1] for x in data]
    else:
        X = [x[:-1] for x in data]
    y = [x[-1] for x in data]

    sm = BorderlineSMOTE(random_state=42, kind="borderline-1")
    X_resample, y_resample = sm.fit_resample(X, y)

    print('Before oversampling, ', 'Pathogenic:', num_pathogenic, 'Benign:', num_benign)
    print('After oversampling, ', 'Pathogenic:', sum(y_resample), 'Benign:', len(y_resample) - sum(y_resample))

    return locationAndId, X, y, X_resample, y_resample

def evaluating(model, true_y, pred_y):
    cm = confusion_matrix(y_true = true_y, y_pred = pred_y).tolist()
    TN = cm[0][0]
    FN = cm[1][0]
    TP = cm[1][1]
    FP = cm[0][1]
    accuracy = round((TP+TN)/(TP+FP+TN+FN), 5)
    precision = round(TP/(TP+FP), 5)
    recall = round(TP/(TP+FN), 5)
    specificity = round(TN/(TN+FP), 5)
    f1_score = round(2*precision*recall/(precision+recall), 5)
    g_mean = round((recall*specificity) ** 0.5, 5)
    mcc = round((TP*TN-FP*FN)/(((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)) ** 0.5), 5)
    return [model, cm, accuracy, precision, recall, specificity, f1_score, g_mean, mcc]

def cross_validation(X, y, X_resample, y_resample, n_splits, test_size):
    GBDT_scores = []
    SVM_scores = []
    LogisticRegression_scores = []
    DecisionTree_scores = []
    KNeighbors_scores = []
    ss = ShuffleSplit(n_splits=n_splits, test_size=test_size, random_state=1)
    k = 0
    for train_index, test_index in ss.split(X_resample):

        train_X = [X_resample[i] for i in train_index]
        train_y = [y_resample[i] for i in train_index]
        test_X = []
        test_y = []

        #sv_id = []
        for i in test_index:
            if X_resample[i] in X:
                test_X.append(X_resample[i])
                test_y.append(y_resample[i])
                #sv_id.append(locationAndId[X.index(X_resample[i])][1])

        # Record the current test set with SV ID for testing SVScore and SVFX performance
        #     test_file = 'test_set_'+str(k)+'.vcf'
        #     with open(test_file, 'w') as write_vcf:
        #         vcf_reader = vcf.Reader(filename=vcf_file)
        #         for record in vcf_reader:
        #             variantID = record.ID
        #             if variantID in sv_id:
        #                 chrom = record.CHROM
        #                 pos = record.POS
        #                 ref = record.REF
        #                 alt = "".join(str(i) for i in record.ALT)
        #                 qual = '.'
        #                 Filter = '.'
        #                 clnsig = record.INFO['CLNSIG']
        #                 clnsig = "".join(str(i) for i in clnsig)
        #                 clnsig = clnsig.split(':')[-1]
        #                 write_vcf.write(str(chrom)+'\t'+str(pos)+'\t'+str(variantID)+'\t'+ref+'\t'+alt+'\t'
        #                                 +qual+'\t'+Filter+'\t'+'CLNSIG='+clnsig[0]+'\n')

        print('     =======================  Round {0}  ======================='.format(k))
        print('The size of the train set is {0}, the number of pathogenic variants is {1}'.format(len(train_X), sum(train_y)))
        print('The size of the test set excludes data from oversampling {0}'.format(len(test_X)))

        tb = pt.PrettyTable()
        tb.field_names = ['Model', 'Confusion Matrix', 'Accuracy', 'Precision', 'Recall', 'Specificity', 'F1-score', 'G-mean', 'MCC']

        model = GradientBoostingClassifier()
        model.fit(train_X, train_y)
        predict_y = model.predict(test_X)
        scores = evaluating('GBDT', test_y, predict_y)
        tb.add_row(scores)
        GBDT_scores.append(scores)

        model = svm.SVC(probability=True)
        model.fit(train_X, train_y)
        predict_y = model.predict(test_X)
        scores = evaluating('SVM', test_y, predict_y)
        tb.add_row(scores)
        SVM_scores.append(scores)

        model = LogisticRegression(max_iter=150)
        model.fit(train_X, train_y)
        predict_y = model.predict(test_X)
        scores = evaluating('Logistic Regression', test_y, predict_y)
        tb.add_row(scores)
        LogisticRegression_scores.append(scores)

        model = DecisionTreeClassifier()
        model.fit(train_X, train_y)
        predict_y = model.predict(test_X)
        scores = evaluating('Decision Tree', test_y, predict_y)
        tb.add_row(scores)
        DecisionTree_scores.append(scores)

        model = KNeighborsClassifier(n_neighbors=3)
        model.fit(train_X, train_y)
        predict_y = model.predict(test_X)
        scores = evaluating('KNeighbors', test_y, predict_y)
        tb.add_row(scores)
        KNeighbors_scores.append(scores)

        print(tb)
        print('\n')
        k += 1
    all_scores = {'GBDT_scores': GBDT_scores, 'SVM_scores': SVM_scores,
                  'LogisticRegression_scores': LogisticRegression_scores,
                  'DecisionTree_scores': DecisionTree_scores, 'KNeighbors_scores': KNeighbors_scores}
    return all_scores

def compute_average(all_scores):
    GBDT_avg = []
    SVM_avg = []
    LogisticRegression_avg = []
    DecisionTree_avg = []
    KNeighbors_avg = []

    GBDT_scores = all_scores['GBDT_scores']
    SVM_scores = all_scores['SVM_scores']
    LogisticRegression_scores = all_scores['LogisticRegression_scores']
    DecisionTree_scores = all_scores['DecisionTree_scores']
    KNeighbors_scores = all_scores['KNeighbors_scores']

    tb = pt.PrettyTable()
    tb.field_names = ['Model', 'Accuracy', 'Precision', 'Recall', 'Specificity', 'F1-score', 'G-mean', 'MCC']

    GBDT_avg.append('GBDT')
    for a in list(np.mean([x[2:] for x in GBDT_scores], axis=0)):
        GBDT_avg.append(round(a, 5))
    tb.add_row(GBDT_avg)

    SVM_avg.append('SVM')
    for a in list(np.mean([x[2:] for x in SVM_scores], axis=0)):
        SVM_avg.append(round(a, 5))
    tb.add_row(SVM_avg)

    LogisticRegression_avg.append('Logistic Regression')
    for a in list(np.mean([x[2:] for x in LogisticRegression_scores], axis=0)):
        LogisticRegression_avg.append(round(a, 5))
    tb.add_row(LogisticRegression_avg)

    DecisionTree_avg.append('Decision Tree')
    for a in list(np.mean([x[2:] for x in DecisionTree_scores], axis=0)):
        DecisionTree_avg.append(round(a, 5))
    tb.add_row(DecisionTree_avg)

    KNeighbors_avg.append('KNeighbors')
    for a in list(np.mean([x[2:] for x in KNeighbors_scores], axis=0)):
        KNeighbors_avg.append(round(a, 5))
    tb.add_row(KNeighbors_avg)

    print(' *********************************** Average ***********************************')
    print(tb)

def train_on_all_data(X_resample, y_resample):
    model = GradientBoostingClassifier()
    model.fit(X_resample, y_resample)
    return model

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Training structural variant pathogenicity prediction model SVPath')
    parser.add_argument('-i', '--input', required=True, help='A CSV file containing variant annotations')
    parser.add_argument('-o', '--output', required=True, help='Path to save the trained model')
    parser.add_argument('-t', '--type', required=True, help='Types of structural variation, DEL, INS and DUP', choices=['DEL', 'INS', 'DUP'])

    args = parser.parse_args()

    features_csv = args.input
    output_path = args.output
    varType = args.type

    locationAndId, X, y, X_resample, y_resample = get_features(features_csv)
    n_splits = 10
    test_size = 0.3

    print('Cross-validating...')
    scores = cross_validation(X, y, X_resample, y_resample, n_splits, test_size)

    print('Calculating the average score...')
    compute_average(scores)

    print('Training the model on all data sets')
    model = train_on_all_data(X_resample, y_resample)
    s = pickle.dumps(model)
    model_name = varType+'.pkl'
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    out_file = os.path.join(output_path, model_name)
    with open(out_file, 'wb+') as f:
        f.write(s)