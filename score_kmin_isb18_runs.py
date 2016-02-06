# PILDA_PROJECT/SCRIPTS/score_isb18_runs.py
# MODIFIED script to score Kmin peptide formatted protein inference runs
# for the ISB18 data sets.
import sys
import argparse
import re
import glob
from collections import defaultdict

"""Maybe unneeded....-JH"""
from functools import partial
from sklearn.metrics import auc

# the ISB18 protein set...

#Modified expression to extract parameters in the format min.{#}.threshold.{#}.ISB18.txt
RUN_PARAMETERS_PATTERN = 'min\.([0-9]+)\.threshold\.([0-9]+)\.ISB18\.txt'
RUN_PARAMETERS_REGEX = re.compile(RUN_PARAMETERS_PATTERN)

def compute_recall(tp, fp, fn):
    denom = tp + fn
    return float(tp)/denom if denom else 0.0
    
def compute_precision(tp, fp, fn):
    denom = tp + fp
    return float(tp)/denom if denom > 0 else 1.0
    
def compute_f1(tp, fp, fn):
    p = compute_precision(tp, fp, fn)
    r = compute_recall(tp, fp, fn)
    denom = p + r
    return 2.0*p*r/denom if denom > 0.0 else 0.0

def compute_roc_auc(xs, ys, positiveYValue):
    """
    Compute area under the receiver operating curve...
    INPUTS:
        xs : predicted floating point values for each sample
        ys : truth for each sample
        positiveYValue : y value to treat as positive
    OUTPUT: computed AUC
    """
    def trapezoid_area(X1, X2, Y1, Y2):
        return 0.50*(X1 - X2)*(Y1 + Y2)

    # arrange the data, sorted by decreasing confidence xs values...
    data = sorted(zip(xs, ys), reverse = True)
    # check for degenerate case, all xs equal...
    if len(xs) == 0 or data[0][0] == data[-1][0]:
        # meaningless, so return 0.50 since no information is gained...
        return 0.50
    # count number of each class...
    P = float(ys.count(positiveYValue))
    N = len(ys) - P
    # check for degenerate cases...
    if P == 0.0 or N == 0.0:
        return 0.0
    # compute ROC points...
    TP = FP = 0.0
    TPprev = FPprev = 0.0
    fprev = None
    A = 0.0
    for x, y in data:
        if x != fprev:
            A += trapezoid_area(FP, FPprev, TP, TPprev)
            FPprev = FP
            TPprev = TP
            fprev = x
        if y == positiveYValue:
            TP += 1.0
        else:
            FP += 1.0
    A += trapezoid_area(N, FPprev, P, TPprev)
    A = A/(P*N)
    return A

def compute_mcc(tp, tn, fp, fn):
    # mcc = Matthews Correlation Coefficient from http://en.wikipedia.org/wiki/Matthews_Correlation_Coefficient
    mccN = float(tp*tn - fp*fn)
    mccD = float((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
    if mccD == 0.0:
        return mccN
    else:
        return mccN/sqrt(mccD)

def compute_average_precision(xs, ys, positiveYValue):
    # this code is rather complicated in order to handle
    # ties with the x value in an unbaised manner. the basic strategy
    # is to use the 'worst case' including all samples of a givnen x value in the precision
    # calculation, and repeating this precision value a number of times equal to the number of
    # positive y samples at that given x value...
    
    # sort data by x value descending...
    data = sorted(zip(xs, ys), key = lambda x : x[0], reverse = True)
    # loop over samples, creating a list of list of samples,
    # where each list of samples have the same x value...
    zss = []
    xstart = None
    for x, y in data:
        # need to start a new list...
        if x != xstart:
            zss.append([])
        # append pair to current list...
        zss[-1].append((x, y))
        xstart = x
    # for each list zs in zss, compute a precision...
    precisions = []
    ny = nx = 0
    for zs in zss:
        dy = 0
        nx += len(zs)
        for (x, y) in zs:
            # count the number of additional positive y
            # samples at this value of x...
            if y == positiveYValue:
                dy += 1
        # if we saw positive y samples at this value of x,
        # compute a precision at this index in the list...
        if dy > 0:
            # compute the average of all the precisions
            # possible when seeing new dy spread across dy..dx additional
            # new samples...
            precision = float(ny + dy) / nx
            precisions.extend([precision]*dy)
        ny += dy    
    # compute the average of the collected precisions...
    if len(precisions) > 0:
        average_precision = sum(precisions)/len(precisions)
    else:
        average_precision = 0.0
    return average_precision

def extract_parameters_from_filepath(filepath, opts):
    if opts.skip_params:
        return dict()

    mo = RUN_PARAMETERS_REGEX.search(filepath)
    if mo:
        return dict(filepath=filepath, min=float(mo.group(1)), threshold=float(mo.group(2)))
    else:
        raise ValueError("Unable to extract parameters from file: %s" % filepath)

def score_predictions_in_filepath(filepath, opts):
    # every prediction belongs to one of these categories:
    # 1. CORRECT ISB18 PROTEIN
    # 2. INCORRECT PROTEIN PREDICTION
    # 3. MISSING ISB18 PROTEIN
    # 4. CONTAMINENT

    runs = defaultdict(lambda : dict())
    with open(filepath, 'rU') as in_file:
        # read the first line header...
        headers = in_file.readline().strip().split('\t')
        # collect the inferences for this file...
        for line in in_file:
            fields = dict(zip(headers, line.strip().split('\t')))
            run_identifier = fields['run_identifier']
            protein = fields['protein']

            #mean_spectral_count = float(fields['mean_spectral_count'])                                    
            """
            mean_spectral_count/pvalue/qvalue not used anywhere else
            changed mean_spectral_count -> total_spectral_count
            """ 
            total_spectral_count = float(fields['total_spectral_count'])                                   
            pvalue = 1#float(fields['pvalue'])                                    
            qvalue = 1#float(fields['qvalue'])

            #if qvalue <= opts.fdr:
                # skip any predictions with too great an fdr...
            runs[run_identifier][protein] = (total_spectral_count, pvalue, qvalue)

        # score the inferences for this file...
        total_correct = total_incorrect = total_missing = total_contaminant = total_undetermined = 0
        for run_identifier in runs:
            correct = incorrect = missing = contaminant = undetermined = other = 0
            for protein in runs[run_identifier]:
                if ("Contaminant" in protein) or protein.startswith("CONT_"):
                    contaminant += 1
                elif "ISB18_Main" in protein:
                    correct += 1
                elif protein == "UNDETERMINED":
                    undetermined += 1
                else:
                    incorrect += 1
            missing = 18 - correct
            total_correct += correct
            total_incorrect += incorrect
            total_missing += missing
            total_contaminant += contaminant
            total_undetermined += undetermined
    return dict(run_count=len(runs), tp=total_correct, fp=total_incorrect, fn=total_missing, total_contaminant=total_contaminant, total_undetermined=total_undetermined)                    

                    
def compute_measures(result, opts):
    tp, fp, fn = result['tp'], result['fp'], result['fn']
    recall = compute_recall(tp, fp, fn)
    precision = compute_precision(tp, fp, fn)
    f1 = compute_f1(tp, fp, fn)
    return dict(recall=recall, precision=precision, f1=f1)


def process_files(filepathlist, opts):

    precision_recalls = defaultdict(lambda : dict())

    parameters_keys = results_keys = measures_keys = None
    for filepath in filepathlist:
        parameters = extract_parameters_from_filepath(filepath, opts)
        results = score_predictions_in_filepath(filepath, opts)
        measures = compute_measures(results, opts)        
        if parameters_keys is None:
            # write header if this is the first file...
            parameters_keys = sorted(parameters)
            results_keys = sorted(results)
            measures_keys = sorted(measures)
            sys.stdout.write("%s\n" % '\t'.join(parameters_keys + results_keys + measures_keys))
        # write out the results...
        values = [str(parameters[k]) for k in parameters_keys]
        values += [str(results[k]) for k in results_keys]
        values += [str(measures[k]) for k in measures_keys]
        sys.stdout.write("%s\n" % '\t'.join(values))


if __name__ == "__main__":
    # configure command line parser...
    parser = argparse.ArgumentParser(description='Command line for %s.' % sys.argv[0])
    parser.add_argument('-in', type=str, help='Path or glob to tab separated text files containing ISB18 protein inference predictions.', required=True)
    parser.add_argument('-fdr', type=float, help='Maximum acceptable false discovery rate (fdr).', required=True)
    parser.add_argument('-skip_params', action='store_true', help="Do not extract parameter settings, just score the proteins.")
    
    opts = parser.parse_args()

    # score all files...
    process_files(glob.glob(getattr(opts, 'in')), opts)
    
