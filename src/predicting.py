import csv
import pandas as pd
import argparse
import pickle

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='SVPath: Used to predict the pathogenicity of deletion structural variations in exons.')
    parser.add_argument('-i', '--input', required=True, help='A csv file containing all feature values.')
    parser.add_argument('-m', '--model', required=True, help='Saved model.')
    parser.add_argument('-o', '--output', required=True, help='A csv file for storing prediction results.')

    args = parser.parse_args()
    input_matrix = args.input
    output_matrix = args.output
    model = args.model

    input_matrix = pd.read_csv(input_matrix)
    cols = input_matrix.columns.to_list()
    input_matrix = pd.DataFrame(input_matrix)
    input_matrix = input_matrix.values.tolist()

    result_matrix = []
    if 'Location' in cols and 'ID' in cols and 'SV_LEN' in cols:
        result_matrix = [x[:3] for x in input_matrix]
        features_matrix = [x[3:] for x in input_matrix]
    else:
        features_matrix = [x for x in input_matrix]

    f = open(model, 'rb')
    model = pickle.load(f)
    f.close()

    predict_result = model.predict(features_matrix)

    num_pathogenic = 0
    num_bengin = 0
    for i in range(len(result_matrix)):
        if predict_result[i] == 1:
            result_matrix[i].append('Pathogenic')
            num_pathogenic += 1
        elif predict_result[i] == 0:
            result_matrix[i].append('Benign')
            num_bengin += 1
        else:
            result_matrix[i].append('unknown')
    print('  ======================================  ')
    print('The number of sv:', len(result_matrix))
    print('Predicting results, pathogenic sv:', num_pathogenic, 'bengin sv:', num_bengin)
    print('  ======================================  ')
    print('Writing prediction results to '+ str(output_matrix) +' file')
    f = open(output_matrix, 'w', encoding='utf-8')
    csv_writer = csv.writer(f)
    csv_head = ['Location', 'ID', 'SV_length', 'SVPath prediction results']
    csv_writer.writerow(csv_head)
    for line in result_matrix:
        try:
            csv_writer.writerow(line)
        except:
            pass
    f.close()