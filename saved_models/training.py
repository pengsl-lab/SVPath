import numpy as np
import pandas as pd
import argparse
import prettytable as pt
from imblearn.over_sampling import BorderlineSMOTE
from sklearn.metrics import confusion_matrix

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
# from sklearn.model_selection import ShuffleSplit
# from sklearn.metrics import matthews_corrcoef, accuracy_score

def get_features(features_csv):
    data = pd.read_csv(features_csv)
    cols = data.columns.to_list()
    data = pd.DataFrame(data)
    data = data.values.tolist()

    print('There are a total of {0} original data'.format(len(data)))
    data = list(set([tuple(t) for t in data]))
    data = [list(v) for v in data]
    print('After removing redundancy, there are {0} remaining'.format(len(data)))

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

    return locationAndId, X, y

def evaluating(model, true_y, pred_y):
    cm = confusion_matrix(y_true = true_y, y_pred = pred_y)
    # TN = cm[0][0]
    # FN = cm[1][0]
    # TP = cm[1][1]
    # FP = cm[0][1]
    TN, FP, FN, TP = cm.ravel()
    accuracy = round((TP+TN)/(TP+FP+TN+FN), 3)
    precision = round(TP/(TP+FP), 3)
    recall = round(TP/(TP+FN), 3)
    specificity = round(TN/(TN+FP), 3)
    f1_score = round(2*precision*recall/(precision+recall), 3)
    g_mean = round((recall*specificity) ** 0.5, 3)
    mcc = round((TP*TN-FP*FN)/(((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)) ** 0.5), 3)
    return [model, cm.tolist(), accuracy, precision, recall, specificity, f1_score, g_mean, mcc]

def cross_validation(locationAndId, X, y, varType):
    GBDT_scores = []
    SVM_scores = []
    LogisticRegression_scores = []
    DecisionTree_scores = []
    KNeighbors_scores = []

    chr_list = []
    train_set_chr = []
    for i in range(1, 25):
        if i == 23:
            chr_list.append('X')
        elif i == 24:
            chr_list.append('Y')
        else:
            chr_list.append((str(i)))

    if varType == 'INS':
        test_set_chr = [['1', '4', '7', '10', '13', '16', '19', '22'],
                        ['2', '5', '8', '11', '14', '17', '20', 'X'],
                        ['3', '6', '9', '12', '15', '18', '21', 'Y']]
    else:
        # test_set_chr = [['1', '9', '17'], ['2', '10', '18'], ['3', '11', '19'], ['4', '12', '20'], ['5', '13', '21'],
        #                 ['6', '14', '22'], ['7', '15', 'X'], ['8', '16', 'Y']]
        test_set_chr = [['1', '7', '13', '19'], ['2', '8', '14', '20'], ['3', '9', '15', '21'],
                        ['4', '10', '16', '22'], ['5', '11', '17', 'X'], ['6', '12', '18', 'Y']]

    for test in test_set_chr:
        train = list(set(chr_list) - set(test))
        train = sorted(train, key=chr_list.index)
        train_set_chr.append(train)

    train_set_X = []
    train_set_y = []
    test_set_X = []
    test_set_y = []

    print('Dividing training set and test set according to chromosomes...')
    for i in range(len(train_set_chr)):
        train_x = []
        train_y = []
        test_x = []
        test_y = []
        for j in range(len(X)):
            if locationAndId[j][0].split(':')[0] in train_set_chr[i]:
                train_x.append(X[j])
                train_y.append(y[j])
            if locationAndId[j][0].split(':')[0] in test_set_chr[i]:
                test_x.append(X[j])
                test_y.append(y[j])
        train_set_X.append(train_x)
        train_set_y.append(train_y)
        test_set_X.append(test_x)
        test_set_y.append(test_y)

    print('The division is completed, divided into {0} groups'.format(len(train_set_X)))
    print('Cross-validating\n')
    for i in range(len(train_set_X)):
        print('     =======================  Round {0}  ======================='.format(i+1))
        print('The current training set consists of the following chromosomes: ', end='')
        for j in range(len(train_set_chr[i])-1):
            print(train_set_chr[i][j] + ', ', end='')
        print(train_set_chr[i][-1])

        train_X = train_set_X[i]
        train_y = train_set_y[i]
        test_X = test_set_X[i]
        test_y = test_set_y[i]

        total = len(train_X)
        num_path = sum(train_y)
        num_benign = total - num_path
        print('Training set sample size: {0}, of which Pathogenic: {1}, Benign:{2}'.format(total, num_path, num_benign))

        # if num_path > 2 * num_benign or  num_benign > 2 * num_path:
        sm = BorderlineSMOTE(random_state=42, kind="borderline-1")
        train_X, train_y = sm.fit_resample(train_X, train_y)

        print('  After oversampling, training set sample size: {0}, of which Pathogenic: {1}, Benign:{2}'.format(len(train_X), sum(train_y), len(train_y) - sum(train_y)))

        print('The current testing set consists of the following chromosomes: ', end='')
        for j in range(len(test_set_chr[i])-1):
            print(test_set_chr[i][j] + ', ', end='')
        print(test_set_chr[i][-1])

        print('The test set sample size is {0}, of which Pathogenic: {1}, Benign: {2} (Original actual variation, '
              'not oversampling)'.format(len(test_X), sum(test_y), len(test_y) - sum(test_y)))

        tb = pt.PrettyTable()
        tb.field_names = ['Model', 'Confusion Matrix', 'Accuracy', 'Precision', 'Recall', 'Specificity', 'F1-score', 'G-mean', 'MCC']

        model = GradientBoostingClassifier()
        model.fit(train_X, train_y)
        predict_y = model.predict(test_X)
        scores = evaluating('GBDT', test_y, predict_y)
        tb.add_row(scores)
        GBDT_scores.append(scores)


        #model = svm.SVC(max_iter=100, probability=True, random_state=1)
        model = svm.NuSVC(probability=True, random_state=1)
        model.fit(train_X, train_y)
        predict_y = model.predict(test_X)
        scores = evaluating('SVM', test_y, predict_y)
        tb.add_row(scores)
        SVM_scores.append(scores)

        model = LogisticRegression(solver='newton-cg', max_iter=100, random_state=1)
        model.fit(train_X, train_y)
        predict_y = model.predict(test_X)
        scores = evaluating('Logistic Regression', test_y, predict_y)
        tb.add_row(scores)
        LogisticRegression_scores.append(scores)

        model = DecisionTreeClassifier(criterion='entropy', splitter='random', max_depth=100)
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
        GBDT_avg.append(round(a, 3))
    tb.add_row(GBDT_avg)

    SVM_avg.append('SVM')
    for a in list(np.mean([x[2:] for x in SVM_scores], axis=0)):
        SVM_avg.append(round(a, 3))
    tb.add_row(SVM_avg)

    LogisticRegression_avg.append('Logistic Regression')
    for a in list(np.mean([x[2:] for x in LogisticRegression_scores], axis=0)):
        LogisticRegression_avg.append(round(a, 3))
    tb.add_row(LogisticRegression_avg)

    DecisionTree_avg.append('Decision Tree')
    for a in list(np.mean([x[2:] for x in DecisionTree_scores], axis=0)):
        DecisionTree_avg.append(round(a, 3))
    tb.add_row(DecisionTree_avg)

    KNeighbors_avg.append('KNeighbors')
    for a in list(np.mean([x[2:] for x in KNeighbors_scores], axis=0)):
        KNeighbors_avg.append(round(a, 3))
    tb.add_row(KNeighbors_avg)

    print(' *********************************** Average ***********************************')
    print(tb)

def train_on_all_data(X, y):
    model = GradientBoostingClassifier()
    sm = BorderlineSMOTE(random_state=42, kind="borderline-1")
    train_X, train_y = sm.fit_resample(X, y)
    model.fit(train_X, train_y)
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

    locationId, X, y = get_features(features_csv)

    print('Cross-validating...')
    scores = cross_validation(locationId, X, y, varType)

    print('Calculating the average score...')
    compute_average(scores)

    print('Training the GBDT model on all data sets')
    model = train_on_all_data(X, y)
    s = pickle.dumps(model)
    model_name = varType+'.pkl'
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    out_file = os.path.join(output_path, model_name)
    with open(out_file, 'wb+') as f:
        f.write(s)
    print('Finished, the pre-trained GBDT model is saved to {0}'.format(os.path.abspath(out_file)))