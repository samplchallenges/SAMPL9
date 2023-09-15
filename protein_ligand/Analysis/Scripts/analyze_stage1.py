#!/usr/bin/env python

# =============================================================================
# GLOBAL IMPORTS
# =============================================================================

import os
import copy
import collections
import pickle

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

from pkganalysis.submission import (SamplSubmission, IgnoredSubmissionError,
                                    load_submissions, plot_correlation)

from pkganalysis.stats import (compute_bootstrap_statistics, calc_confusion_matrix, balanced_accuracy,
                              accuracy, sensitivity, specificity, precision,
                            TP, TN, FP, FN, cohen_kappa, mcc, j)

from pkganalysis.stats import (cohen_kappa, mcc, j)

#import sklearn metrics and other modules for random predictions
from sklearn.metrics import confusion_matrix, balanced_accuracy_score, f1_score, recall_score, precision_score, precision_recall_curve, roc_curve, cohen_kappa_score, matthews_corrcoef

#for plotting confusion matrix
from sklearn.metrics import ConfusionMatrixDisplay

from matplotlib.colors import LogNorm
from matplotlib.ticker import MaxNLocator

from statistics import mean
import random


# =============================================================================
# TO DO NOTES
# =============================================================================

# Script adapted from SAMPL6/SAMPL7 (Andrea Rizzi and David Mobley), and SAMPL7 Protein-Ligand (Harold Grosjean). 

# UPDATES

# - (DONE) Check submission files, if method name of submission is not unique, assign one
# - (DONE) Udate to include "ranked only" and "Non ranked", and "All" (ranked and non-ranked) analysis output separate
# - (DONE) Include: Overal percent of true actives, false active, true non-actives, and false non-actives
# - (DONE) Include: Bindary statistics: See SAMPL7 Protein-ligand Analysis
# - (DONE) Make a new experimental results hits_verification.csv file. Need to include True values to the inhibitors, include non-inhibitors
# and include False values for the non-inhibitors, then arrange these based on name. 
# - (DONE) For submissions which only submitted 'inhibitors' and their values as True, fill in all other SAMPL9 dataset
# ligand names and give them a value of False (non-inhibitors).
# - (DONE) Add Youdens J Index, ROC, MCC (sklearn) 
# - (DONE) Generate Confusion matrix plots for each submission, 
# - (DONE) Generate plots for each metric calculated
# - (DONE) Generate some "Random Predictions" with known p_active = 2083/45543
# - (DONE) Include random performance to each of the metric plots


# TO DO / WIP: Update analyze_stage1.py 
# - (DONE) Add MCC/Phi, and cohens kappa metrics to stats.py
# - Recall vs Precision curve to stats.py

# - (DONE) Include the generation of random predictions with p_active in plots of metrics. 
# - (DONE) Generate Confusion matrix plots for each submission (also for each ranking type (i.e Ranked, Non-Ranked, Combined))

# - (DONE) Generate plots for each metric calculated with random compmarison, and save to respective directory

# - (DONE) Change confusion matrices to log scale.


# =============================================================================
# CONSTANTS
# =============================================================================

# Paths to input data.
STAGE_1_SUBMISSIONS_DIR_PATH = '/Users/amezcum1/Desktop/SAMPL9/protein_ligand/Analysis/Submissions/'
EXPERIMENTAL_DATA_FILE_PATH = '/Users/amezcum1/Desktop/SAMPL9/protein_ligand/Analysis/Scripts/hits_verification.csv'
USER_MAP_FILE_PATH = '/Users/amezcum1/Desktop/SAMPL9/protein_ligand/Analysis/SAMPL9-user-map.csv'
FRAGMENTS_FILE_PATH = '/Users/amezcum1/Desktop/SAMPL9/protein_ligand/NCATS_experimental_data/nanoluc_compounds_tranche1.csv'

# =============================================================================
# MAIN CHALLENGE PROTEIN-LIGAND SUBMISSION
# =============================================================================

class Stage1Submission(SamplSubmission):
    """A submission for the main protein-ligand challenge.

    Parameters
    ----------
    file_path : str
        The path to the submission file.

    Raises
    ------
    IgnoredSubmission
        If the submission ID is among the ignored submissions.

    """

    # The IDs of the submissions used for testing the validation. Should be strings of submission IDs

    TEST_SUBMISSION_SIDS = {}
    # The IDs of submissions for reference calculations. Should be strings of submission IDs
    REF_SUBMISSION_SIDS = [3, 4]  

    # Section of the submission file.
    SECTIONS = {'Predictions', 'Participant name', 'Participant organization', 'Name', 'Software', 'Method', 'Category', 'Ranked'}


    # Sections in CSV format with kwargs to pass to pandas.read_csv().

    CSV_SECTIONS = {
        'Predictions': {'names': ('Fragment ID', 'Site 1'),
                        'index_col': 'Fragment ID'}
    }

    RENAME_METHODS = {}

    def __init__(self, file_path, user_map):
        super().__init__(file_path, user_map)

        file_name = os.path.splitext(os.path.basename(file_path))[0]
        self.file_name = file_name

        # Load predictions.
        sections = self._load_sections(file_path)  # From parent-class.
        self.data = sections['Predictions']  # This is a list
        self.data = pd.DataFrame(data=self.data) # Now a DataFrame
        try:
            self.name = self.RENAME_METHODS[sections['Name'][0]]
        except KeyError:
            self.name = sections['Name'][0]


        # Store participant name, organization, method category
        self.participant = sections['Participant name'][0].strip()
        self.category = sections['Category'][0].strip()
        self.organization = sections['Participant organization'][0].strip()
        self.ranked = sections['Ranked'][0].strip() =='True'

        # Check if this is a test submission.
        if self.sid in self.TEST_SUBMISSION_SIDS:
            raise IgnoredSubmissionError('This submission has been used for tests.')

        # Check if this is a reference submission
        self.reference_submission = False
        if self.sid in self.REF_SUBMISSION_SIDS:
            self.reference_submission = True

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

class Stage1SubmissionCollection:
    """A collection of Stage 1 submissions."""


    def __init__(self, submissions, experimental_data, output_directory_path, stage1_submission_collection_file_path, ignore_refcalcs = False):


        # Check if submission collection file already exists.
        if os.path.isfile(stage1_submission_collection_file_path):
            print("Analysis will be done using the existing submission collection file: {}".format(stage1_submission_collection_file_path))

            self.data = pd.read_csv(stage1_submission_collection_file_path)
            print("\n SubmissionCollection: \n")
            print(self.data)

            # Populate submission.data dataframes parsing sections of collection file.
            for submission in submissions:

                # To ignore reference calculations, when necessary
                if submission.reference_submission and ignore_refcalcs:
                    continue

                df_collection_of_each_submission = self.data.loc[self.data["SID"] == int(submission.sid)]

                # Transform into Pandas DataFrame.
                submission.data = pd.DataFrame()
                submission.data["Site 1"] = df_collection_of_each_submission["Site 1 (pred)"]
                submission.data["Fragment ID"] = df_collection_of_each_submission["Fragment ID"]

                submission.data.set_index("Fragment ID", inplace=True)

            # Transform into Pandas DataFrame.
            self.output_directory_path = output_directory_path


        else: # Build collection dataframe from the beginning.
            # Build full stage 1 collection table.
            data = []

            # Submissions for stage 1.
            for submission in submissions:
                if submission.reference_submission and ignore_refcalcs:
                    continue

                #print("submission.sid:\n", submission.sid)
                #print("submission.name:\n", submission.name)
                #print("submission.data:\n", submission.data)


                for fragment_ID, series in submission.data.iterrows():
                    #print("fragment_ID:", fragment_ID)
                    #print("series:\n", series)
                    #site1_pred = series["Site 1"]
                    #print("site1_pred: ", site1_pred)

                    # Predicted data
                    site1_pred = series["Site 1"]

                    # Experimental data
                    site1_exp = experimental_data.loc[fragment_ID, 'Site 1']

                    data.append({
                        'SID': submission.sid,  # Previously receipt_ID
                        'Participant': submission.participant,
                        'Organization': submission.organization,
                        'Name': submission.name,
                        'Category': submission.category,
                        'Ranked': submission.ranked,
                        'Fragment ID': fragment_ID,
                        'Site 1 (pred)': site1_pred,
                        'Site 1 (exp)': site1_exp,
                    })

            # Transform into Pandas DataFrame.
            self.data = pd.DataFrame(data=data)
            self.output_directory_path = output_directory_path

            print("\n SubmissionCollection: \n")
            print(self.data)

            # Create general output directory.
            os.makedirs(self.output_directory_path, exist_ok=True)

            # Save collection.data dataframe in a CSV file.
            self.data.to_csv(stage1_submission_collection_file_path, index=False)

            print("Stage1 submission collection file generated:\n", stage1_submission_collection_file_path)
       
    #Function to complete submission predictions (i.e., if participant only includes inhibitors, we add missing molecules as non-inhibitors)
    def complete_predictions_with_missing_fragments(self, fragments_file_path, submission_collection_file_path,
                                                    ranking):
        """ Adds missing non-binder predictions to collection. Assumes
        Parameters
        ----------
        fragments_file_path: Path to CSV file with fragment IDs in the first column
        filter_nonranked: filters out non ranked results
        """

        # Read fragments file to extract full list of screened fragment list IDs
        # skip the first row for here since these are "headers"
        #(i.e, Sampl_ID, Site 1, and SMILES end up as first row)
        fragments_data = pd.read_csv(fragments_file_path, names=["Fragment ID", "Site 1", "SMILES"], skiprows=[0])
        fragment_IDs = fragments_data["Fragment ID"].values

        # Rebuild full stage 1 collection table, adding missing fragments (assumed they are non-inhibitors i.e. False).
        data = []

        # Submissions for stage 1.
        for submission in submissions:
            if submission.reference_submission and ignore_refcalcs:
                continue

            # print("submission.sid:\n", submission.sid)
            # print("submission.name:\n", submission.name)
            # print("submission.data:\n", submission.data)

            for screened_fragment_ID in fragment_IDs:
                #print("screened_fragment_ID:", screened_fragment_ID)

                # Check if screened fragment is in submission
                submitted_fragment_IDs = set(submission.data.index.values)
               # print("submitted_fragment_IDs:\n", submitted_fragment_IDs)

                # If screened fragment ID is already in the submission set, take prediction records from the submission
                if screened_fragment_ID in submitted_fragment_IDs:
                    #print("Already submitted.")

                    series = submission.data.loc[screened_fragment_ID,:]

                    # Predicted data
                    site1_pred = series["Site 1"]

                # If screened fragment ID was not in the submitted prediction set, add a non-binder prediction to data
                else:
                    #print("Not submitted.")

                    # Predicted data
                    site1_pred = "False"

                # Experimental data
                site1_exp = experimental_data.loc[screened_fragment_ID, 'Site 1']

                data.append({
                   'SID': submission.sid,  # Previously receipt_ID
                    'Participant': submission.participant,
                    'Organization': submission.organization,
                    'Name': submission.name,
                    'Category': submission.category,
                    'Ranked': submission.ranked,
                    'Fragment ID': screened_fragment_ID,
                    'Site 1 (pred)': site1_pred,
                    'Site 1 (exp)': site1_exp
                    })

            #print("data:\n", data)

            # Transform into Pandas DataFrame.
            self.data = pd.DataFrame(data=data)

            #filter ranked, non-ranked, and ranked and non-ranked submission analysis into separate directories
            if ranking == ('Ranked_and_non-ranked' or None):
                self.data = self.data
            if ranking == 'Ranked':
                self.data = self.data[self.data.Ranked == True]
            if ranking == 'Non-ranked':
                self.data = self.data[self.data.Ranked == False]

            print("\n SubmissionCollection: \n")
            print(self.data)

            # Create general output directory.
            os.makedirs(self.output_directory_path, exist_ok=True)

            # Save completed collection.data dataframe in a CSV file.
            self.data.to_csv(submission_collection_file_path, index=False)
            print("Stage1 submission collection file updated with missing predictions:\n", submission_collection_file_path)



    # TO DO: The following function does not pertain to this challenge/needs updating if we even retain
    @staticmethod
    def _assign_paper_method_name(name):
        return name

        
    #generate statistics tables 
    def generate_statistics_tables(self, stats_funcs, subdirectory_path, groupby, site,
                                   extra_fields=None, sort_stat=None,
                                   ordering_functions=None, latex_header_conversions=None,
                                   caption=''):
        """Generate statistics tables in CSV, JSON, and LaTex format.
        Parameters
        ----------
        groupby : str
            The name of the data column to be used to compute the statistics.
            For example, 'name' to obtain statistics about individual methods,
            'system_id' to compute statistics by molecules.
        ordering_functions : dict
            Dictionary statistic_name -> ordering_function(stats), where
            ordering_function determines how to rank the the groups by
            statistics.
        """

        if extra_fields is None:
            extra_fields = []

        def escape(s):
        #    return s.replace('_', '\_')
            return s

        extra_fields_latex = [escape(extra_field) for extra_field in extra_fields]

        file_base_name = 'statistics'
        directory_path = os.path.join(self.output_directory_path, subdirectory_path)

        stats_names, stats_funcs = zip(*stats_funcs.items())
        ci_suffixes = ('', '_lower_bound', '_upper_bound')

        # Compute or read the bootstrap statistics from the cache.
        cache_file_path = os.path.join(self.output_directory_path, 'bootstrap_distributions.p')
        all_bootstrap_statistics = self._get_bootstrap_statistics(groupby, site, stats_names, stats_funcs,
                                                                  cache_file_path=cache_file_path)
        # Collect the records for the DataFrames.
        statistics_csv = []
        statistics_latex = []

        groups = self.data[groupby].unique()
        for i, group in enumerate(groups):
            print('\rGenerating bootstrap statistics tables for {} {} ({}/{})'
                  ''.format(groupby, group, i+1, len(groups)), end='')

            # Isolate bootstrap statistics.
            bootstrap_statistics = all_bootstrap_statistics[group]

            # Select the group.
            data = self.data[self.data[groupby] == group]


            # Isolate the extra field.
            group_fields = {}
            latex_group_fields = {}
            for extra_field, extra_field_latex in zip(extra_fields, extra_fields_latex):
                assert len(data[extra_field].unique()) == 1
                extra_field_value = data[extra_field].values[0]
                group_fields[extra_field] = extra_field_value
                latex_group_fields[extra_field_latex] = escape(extra_field_value)

            record_csv = {}
            record_latex = {}
            for stats_name, (stats, (lower_bound, upper_bound), bootstrap_samples) in bootstrap_statistics.items():
                # For CSV and JSON we put confidence interval in separate columns.
                for suffix, info in zip(ci_suffixes, [stats, lower_bound, upper_bound]):
                    record_csv[stats_name + suffix] = info

                # For the PDF, print bootstrap CI in the same column.
                stats_name_latex = latex_header_conversions.get(stats_name, stats_name)
                record_latex[stats_name_latex] = '{:.2f} [{:.2f}, {:.2f}]'.format(stats, lower_bound, upper_bound)

            statistics_csv.append({'ID': group, **group_fields, **record_csv})
            statistics_latex.append({'ID': escape(group), **latex_group_fields,
                                     **record_latex})

        # Convert dictionary to Dataframe to create tables/plots easily.
        statistics_csv = pd.DataFrame(statistics_csv)
        statistics_csv.set_index('ID', inplace=True)
        statistics_latex = pd.DataFrame(statistics_latex)

        # Sort by the given statistics.
        if sort_stat is not None:
            ordering_function = ordering_functions.get(sort_stat, lambda x: x)
            order = sorted(statistics_csv[sort_stat].items(), key=lambda x: ordering_function(x[1]))
            order = [k for k, value in order]
            statistics_csv = statistics_csv.reindex(order)
            latex_order = [escape(k) for k in order]
            statistics_latex.ID = statistics_latex.ID.astype('category')
            statistics_latex.ID.cat.set_categories(latex_order, inplace=True)
            statistics_latex.sort_values(by='ID', inplace=True)

        # Reorder columns that were scrambled by going through dictionaries.
        stats_names_csv = [name + suffix for name in stats_names for suffix in ci_suffixes]
        stats_names_latex = [latex_header_conversions.get(name, name) for name in stats_names]
        statistics_csv = statistics_csv[extra_fields + stats_names_csv]
        statistics_latex = statistics_latex[['ID'] + extra_fields_latex + stats_names_latex]

        # Create CSV and JSON tables (correct LaTex syntax in column names).
        os.makedirs(directory_path, exist_ok=True)
        file_base_path = os.path.join(directory_path, file_base_name)
        with open(file_base_path + '.csv', 'w') as f:
            statistics_csv.to_csv(f)
        with open(file_base_path + '.json', 'w') as f:
            statistics_csv.to_json(f, orient='index')

        # Create LaTex table.
        latex_directory_path = os.path.join(directory_path, file_base_name + 'LaTex')
        os.makedirs(latex_directory_path, exist_ok=True)
        with open(os.path.join(latex_directory_path, file_base_name + '.tex'), 'w') as f:
            f.write('\\documentclass[8pt]{article}\n'
                    '\\usepackage[a4paper,margin=0.2in,tmargin=0.5in,bmargin=0.5in,landscape]{geometry}\n'
                    '\\usepackage{booktabs}\n'
                    '\\usepackage{longtable}\n'
                    '\\pagenumbering{gobble}\n'
                    '\\begin{document}\n'
                    '\\begin{center}\n'
                    '\\begin{footnotesize}\n')
            statistics_latex.to_latex(f, column_format='|' + 'c'*(2 + len(stats_funcs)) + '|',
                                      escape=False, index=False, longtable=True, bold_rows=True)
            f.write('\end{footnotesize}\n'
                    '\end{center}\n')
            f.write(caption + '\n')
            f.write('\end{document}\n')

        
        #generate some random predictions and compute metrics
        #define predictions to make (k), and number of samples (iterations)
        k=45543
        iterations=1000
        
        #redefine experimental file path for this function
        EXPERIMENTAL_FILE_PATH='/Users/amezcum1/Desktop/SAMPL9/protein_ligand/Analysis/Scripts/hits_verification.csv'
        
        print("Generating {} Samples of Random Predictions".format(iterations))
        
        p_active = (2083/45543) *100
        p_inactive = 100 - p_active
        
        #generate 1000 samples of random predictions
        random_predictions = []
        
        for i in range(iterations):
            random_predictions.append(random.sample(random.choices(range(2), weights=(p_inactive, p_active), k=k), k=k))
            
        #make temporary df of predictions
        random_df = pd.DataFrame(random_predictions)
        #invert rows and columns
        random_df = random_df.transpose()
        
        #add experimental data to random predictions df
        exp_df = pd.read_csv(EXPERIMENTAL_FILE_PATH, 
                            names=['Fragment ID', 'Site 1 (exp)', 'SMILES'])
        
        #drop first row, and drop smiles column of exp df
        exp_df = exp_df.drop(0, axis=0)
        exp_df = exp_df.drop('SMILES', axis=1)
        
        #reset index of exp_df to match random df
        exp_df = exp_df.reset_index(drop=True) #drop=True avoids adding old index as column)
        
        #extract the fragment id and exp value columns from exp_df, then add to random_df
        names = exp_df['Fragment ID']
        exp_values = exp_df['Site 1 (exp)']
        
        random_df = random_df.join(names)
        random_df = random_df.join(exp_values)
        
        #convert all random predictions to booleans
        for column in random_df.columns[:-2]: #all columns except last two (names and exp values)
            random_df[column] = random_df[column].astype(bool)
        
        #use mapping to convert exp column from strings to booleans. Otherwise can't compute metrics
        d = {'True': True, 'False': False}
        random_df['Site 1 (exp)'] = random_df['Site 1 (exp)'].map(d)
        
        #define some metrics to compute
        def specificity_score(tn, tp):
            specificity = tn / (tn + tp)
            return specificity
        
        def j(sensitivity, specificity):
            j = sensitivity + specificity - 1
            return j
        
        #compute metrics
        print("Computing metrics for random predictions")
        random_sensitivities = []
#         random_f1s = []
        random_precisions = []
        random_specificities = []
        random_balanced_accuracies = []
        random_js = []
        random_mccs = []
        random_ks = []
        
        for column in random_df.columns[:-2]:
            #compute confusion matrices
            tn, fp, fn, tp = confusion_matrix(random_df['Site 1 (exp)'], random_df[column], labels=[False, True]).ravel()
            
            #compute each metric for each iteration of random predictions
#             f1 = f1_score(random_df['Site 1 (exp)'], df[column], average='binary')
            sensitivity = recall_score(random_df['Site 1 (exp)'], random_df[column])
            precision = precision_score(random_df['Site 1 (exp)'], random_df[column])
            specificity = specificity_score(tn, fp)
            balanced_accuracy = balanced_accuracy_score(random_df['Site 1 (exp)'], random_df[column])
            js = j(sensitivity, specificity)
            mcc = matthews_corrcoef(random_df['Site 1 (exp)'], random_df[column])
            k = cohen_kappa_score(random_df['Site 1 (exp)'], random_df[column])

            random_sensitivities.append(sensitivity)
#             random_f1s.append(f1)
            random_precisions.append(precision)
            random_specificities.append(specificity)
            random_balanced_accuracies.append(balanced_accuracy)
            random_js.append(js)
            random_mccs.append(mcc)
            random_ks.append(k)
        
        #compute average for each metric
        average_random_sensitivity = mean(random_sensitivities)
#         average_random_f1 = mean(random_f1s)
        average_random_precision = mean(random_precisions)
        average_random_specificity = mean(random_specificities)
        average_random_balanced_accuracy = mean(random_balanced_accuracies)
        average_random_j = mean(random_js)
        average_random_mcc = mean(random_mccs)
        average_random_k = mean(random_ks)
        print("Done computing metrics for random predictions!")
        
        
        #load relevant stats into dataframe
        statistics_df = pd.read_csv(file_base_path + '.csv')
        

        #highlight the null model red in the plots if its there. 
        bar_color = [('red' if b==6 else 'C0') for b in statistics_df['ID']]
        
        
        #plot recall (Sensitivity or True positive rate (TPR))
        print("Generating metric plots for {}".format(ranking))
        plt.figure(figsize=(30, 20))
        statistics_df.plot(x='ID', y='Sensitivity', kind='barh',
                label='_', #removes label
                xerr=[(statistics_df['Sensitivity'] - statistics_df['Sensitivity_lower_bound']),
                      (statistics_df['Sensitivity_upper_bound'] - statistics_df['Sensitivity'])],
                color=bar_color)
        plt.axvline(x=average_random_sensitivity, color='orange', linestyle='--', label='Random') #random predictions
        plt.xlim(0,0.30) #made limit smaller, so we can see differences between methods
        plt.ylabel("SID", fontsize=16)
        plt.xlabel("Rate (TP / (TP+FN)",fontsize=16)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.title("Stage1 - Recall (TPR) for {}".format(ranking), loc='center', fontsize=20)
        plt.legend(bbox_to_anchor=(1.0, 1), loc='upper left')
        plt.savefig('{}/Recall_of_{}_predictions.png'.format(OUTPUT_DIRECTORY_PATH_SPECIFIC, ranking), bbox_inches = 'tight')
        plt.close()

        # plot specificity (true negative rate (TNR)
        plt.figure(figsize=(30,20))
        statistics_df.plot(x='ID', y='Specificity', kind='barh',
                label='_',
                xerr=[(statistics_df['Specificity'] - statistics_df['Specificity_lower_bound']),
                      (statistics_df['Specificity_upper_bound'] - statistics_df['Specificity'])],
                color=bar_color)
        plt.axvline(x=average_random_specificity, color='orange', linestyle='--', label='Random')
        plt.xlim(0, 1)
        plt.ylabel("SID", fontsize=16)
        plt.xlabel("Rate (TN / (TN+FP)", fontsize=16)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.title('Stage1 - Specificity (TNR) for {}'.format(ranking), loc='center', fontsize=20)
        plt.legend(bbox_to_anchor=(1.0, 1), loc='upper left')
        plt.savefig('{}/Specificity_of_{}_predictions.png'.format(OUTPUT_DIRECTORY_PATH_SPECIFIC, ranking), bbox_inches = 'tight')
        plt.close()

        # Plot precision (correctly classified positives (TP / (TP+FP))
        plt.figure(figsize=(30,20))
        statistics_df.plot(x='ID', y='Precision', kind='barh',
                label='_',
                xerr=[(statistics_df['Precision'] - statistics_df['Precision_lower_bound']),
                      (statistics_df['Precision_upper_bound'] - statistics_df['Precision'])],
                color=bar_color)
        plt.axvline(x=average_random_precision, color='orange', linestyle='--', label='Random')
        plt.xlim(0, 0.3) # made limit smaller to see differences between methods
        plt.ylabel("SID", fontsize=16)
        plt.xlabel("Rate (TP / (TP+FP)", fontsize=16)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.title('Stage1 - Precision for {}'.format(ranking), loc='center', fontsize=20)
        plt.legend(bbox_to_anchor=(1.0, 1), loc='upper left')
        plt.savefig('{}/Precision_of_{}_predictions.png'.format(OUTPUT_DIRECTORY_PATH_SPECIFIC, ranking), bbox_inches = 'tight')
        plt.close()
        
        #plot F1 Score
#         statistics_df.plot(x='ID', y='F1 Score', kind='barh',
#                 label='_',
#                 xerr=[(statistics_df['F1 Score'] - statistics_df['F1 Score_lower_bound']),
#                       (statistics_df['F1 Score_upper_bound'] - statistics_df['F1 Score'])],
#                 color=bar_color)
#         plt.axvline(x=average_random_precision, color='orange', linestyle='--', label='Random')
#         plt.xlim(0, 0.3) # made limit smaller to see differences between methods
#         plt.ylabel("SID", fontsize=20)
#         plt.xlabel("F1", fontsize=20)
#         plt.xticks(fontsize=15)
#         plt.yticks(fontsize=15)
#         plt.title('Stage1 - Harmonic Mean of Precision and Recall for {}'.format(ranking), loc='center', fontsize=25)
#         plt.legend(bbox_to_anchor=(1.0, 1), loc='upper left')
#         plt.savefig('{}/F1_of_{}_predictions.png'.format(OUTPUT_DIRECTORY_PATH_SPECIFIC, ranking))
#         plt.close()
        
        #plot youdens j (balanced accuracy)
        plt.figure(figsize=(30,20))
        statistics_df.plot(x='ID', y='Youdens J', kind='barh',
                label='_',
                xerr=[(statistics_df['Youdens J'] - statistics_df['Youdens J_lower_bound']),
                      (statistics_df['Youdens J_upper_bound'] - statistics_df['Youdens J'])],
                color=bar_color)
        plt.axvline(x=average_random_j, color='orange', linestyle='--', label='Random')
        plt.xlim(-0.05, 0.2) # made limit smaller to see differences between methods
        plt.ylabel("SID", fontsize=16)
        plt.xlabel("J", fontsize=16)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.title('Stage1 - Youdens J for {}'.format(ranking), loc='center', fontsize=20)
        plt.legend(bbox_to_anchor=(1.0, 1), loc='upper left')
        plt.savefig('{}/Youdens_of_{}_predictions.png'.format(OUTPUT_DIRECTORY_PATH_SPECIFIC, ranking), bbox_inches = 'tight')
        plt.close()
        
        #Plot MCC (mathhews correlation coefficient)
        plt.figure(figsize=(30,20))
        statistics_df.plot(x='ID', y='MCC', kind='barh',
                label='_',
                xerr=[(statistics_df['MCC'] - statistics_df['MCC_lower_bound']),
                      (statistics_df['MCC_upper_bound'] - statistics_df['MCC'])],
                color=bar_color)
        plt.axvline(x=average_random_mcc, color='orange', linestyle='--', label='Random')
        plt.xlim(-0.05, 0.2) # made limit smaller to see differences between methods
        plt.ylabel("SID", fontsize=16)
        plt.xlabel("MCC", fontsize=16)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.title('Stage1 - Matthews Correlation Coefficient for {}'.format(ranking), loc='center', fontsize=20)
        plt.legend(bbox_to_anchor=(1.0, 1), loc='upper left')
        plt.savefig('{}/MCC_of_{}_predictions.png'.format(OUTPUT_DIRECTORY_PATH_SPECIFIC, ranking), bbox_inches = 'tight')
        plt.close()
        
        #Plot Cohen Kappa ()
        plt.figure(figsize=(30,20))
        statistics_df.plot(x='ID', y='Cohen Kappa', kind='barh',
                label='_',
                xerr=[(statistics_df['Cohen Kappa'] - statistics_df['Cohen Kappa_lower_bound']),
                      (statistics_df['Cohen Kappa_upper_bound'] - statistics_df['Cohen Kappa'])],
                color=bar_color)
        plt.axvline(x=average_random_k, color='orange', linestyle='--', label='Random')
        plt.xlim(-0.05, 0.2) # made limit smaller to see differences between methods
        plt.ylabel("SID", fontsize=16)
        plt.xlabel("Cohen Kappa Score", fontsize=16)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.title('Stage1 - Cohen Kappa for {}'.format(ranking), loc='center', fontsize=20)
        plt.legend(bbox_to_anchor=(1.0, 1), loc='upper left')
        plt.savefig('{}/Cohen_Kappa_of_{}_predictions.png'.format(OUTPUT_DIRECTORY_PATH_SPECIFIC, ranking), bbox_inches = 'tight')
        plt.close()
        
        #plot balanced accuracies (should be similar or same as j)
        plt.figure(figsize=(30,20))
        statistics_df.plot(x='ID', y='Balanced Accuracy', kind='barh',
                label='_',
                xerr=[(statistics_df['Balanced Accuracy'] - statistics_df['Balanced Accuracy_lower_bound']),
                      (statistics_df['Balanced Accuracy_upper_bound'] - statistics_df['Balanced Accuracy'])],
                color=bar_color)
        plt.axvline(x=average_random_balanced_accuracy, color='orange', linestyle='--', label='Random')
        plt.xlim(0, 1) 
        plt.ylabel("SID", fontsize=16)
        plt.xlabel("Rate ", fontsize=16)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.title('Stage1 - Balanced Accuracy for {}'.format(ranking), loc='center', fontsize=20)
        plt.legend(bbox_to_anchor=(1.0, 1), loc='upper left')
        plt.savefig('{}/Balanced_Accuracy_of_{}_predictions.png'.format(OUTPUT_DIRECTORY_PATH_SPECIFIC, ranking), bbox_inches = 'tight')
        plt.close()
        print("Finished generating metric plots for {}".format(ranking))



    def _get_bootstrap_statistics(self, groupby, site, stats_names, stats_funcs, cache_file_path):
        """Generate the bootstrap distributions of all groups and cache them.
        If cached values are found on disk, the distributions are not recomputed.
        Returns
        -------
        all_bootstrap_statistics : collections.OrderedDict
            group -> {stats_name -> (statistics, confidence_interval, bootstrap_samples)}
            confidence_interval is a pair (lower_bound, upper_bound), and bootstrap_samples
            are the (ordered) bootstrap statistics used to compute the confidence interval.
        """


        # Identify all the groups (e.g. methods/molecules).
        groups = self.data[groupby].unique()

        # Initialize returned value. The OrderedDict maintains the order of statistics.
        all_bootstrap_statistics = collections.OrderedDict([(name, None) for name in stats_names])
        all_bootstrap_statistics = collections.OrderedDict(
            [(group, copy.deepcopy(all_bootstrap_statistics)) for group in groups]
        )

        # Load the statistics that we have already computed.
        try:
            with open(cache_file_path, 'rb') as f:
                print('Loading cached bootstrap distributions from {}'.format(cache_file_path))
                cached_bootstrap_statistics = pickle.load(f)
        except FileNotFoundError:
            cached_bootstrap_statistics = None

        # Create a map from paper method name to submission method name.
        try:
            paper_to_submission_name = {self._assign_paper_method_name(submission_name): submission_name
                                        for submission_name in cached_bootstrap_statistics}
        except (UnboundLocalError, TypeError):
            # cached_bootstrap_statistics is None or group is not a method.
            paper_to_submission_name = {}

        cache_updated = False
        for i, (group, group_bootstrap_statistics) in enumerate(all_bootstrap_statistics.items()):
            # Check which statistics we still need to compute for this group.
            if cached_bootstrap_statistics is not None:
                group_stats_names = []
                group_stats_funcs = []
                for stats_name, stats_func in zip(stats_names, stats_funcs):
                    try:
                        all_bootstrap_statistics[group][stats_name] = cached_bootstrap_statistics[group][stats_name]
                    except KeyError:
                        try:
                            # method_name = self._assign_paper_method_name(group)
                            method_name = paper_to_submission_name[group]
                            all_bootstrap_statistics[group][stats_name] = cached_bootstrap_statistics[method_name][
                                stats_name]
                        except KeyError:
                            group_stats_names.append(stats_name)
                            group_stats_funcs.append(stats_func)
            else:
                # Compute everything.
                group_stats_names = stats_names
                group_stats_funcs = stats_funcs

            if len(group_stats_names) == 0:
                continue
            cache_updated = True  # Update the cache on disk later.

            print('\rGenerating bootstrap statistics for {} {} ({}/{})'
                  ''.format(groupby, group, i + 1, len(groups)), end='')

            # Select the group data.
            data = self.data[self.data[groupby] == group]
            print(type(data))
            #print("data:\n", data)

            # Compute bootstrap statistics.
            # Modify here to do per-site statistic
            #loop iterates over: sites = ["site-1", "site-2", "site-3", "site-4", "all-sites"]

            if site == "site-1":
                data = data[['Site 1 (exp)', 'Site 1 (pred)']]

            #mixutre of strings and bools present in df. Converts everything into bools.
            data = data.replace([0.00000, 0, 'False', 'FALSE', 'false'], False)
            data = data.replace([1.00000, 1, 'True', 'TRUE', 'true'], True)

            new_bootstrap_statistics = compute_bootstrap_statistics(data.to_numpy(), group_stats_funcs, sems=None,
                                                                    n_bootstrap_samples=10000) #10000

            #new_bootstrap_statistics returns correct values
            # Update the returned value with the statistics just computed.
            new_boostrap_statistics = {group_stats_names[i]: new_bootstrap_statistics[i]
                                       for i in range(len(group_stats_funcs))}


            print(new_bootstrap_statistics)
            group_bootstrap_statistics.update(new_boostrap_statistics)


        # Cache the computed statistics on disk. Create output directory if necessary.
        if cache_updated:
            os.makedirs(os.path.dirname(cache_file_path), exist_ok=True)
            with open(cache_file_path, 'wb') as f:
                pickle.dump(all_bootstrap_statistics, f)

        return all_bootstrap_statistics


# =============================================================================
# MAIN
# =============================================================================

if __name__ == '__main__':

    # Read experimental data.
    with open(EXPERIMENTAL_DATA_FILE_PATH, 'r') as f:
        experimental_data = pd.read_csv(f, sep=',', index_col='Sampl_ID', header=0)
        print("Experimental data of stage 1:\n",experimental_data)


    #Import user map. USER_MAP CAN'T HAVE ANY SPACE BETWEEN "TABS". 
    # example: 1,file_name. (sid, filename of submission)
    try:
        with open(USER_MAP_FILE_PATH, 'r') as f:
            user_map = pd.read_csv(f)
            print("user_map:\n", user_map)
    except FileNotFoundError:
        user_map=None
        print("Warning: No user map found.")

    # Configuration: statistics to compute.
    # ADD BINARY STATISTICS TO COMPUTE. (I.E. ROC-AUC, MCC/Phi, cohen kappa etc., ? )
    stats_funcs = collections.OrderedDict([
        ('Sensitivity', sensitivity),
        ('Specificity', specificity),
        ('Precision', precision),
        ('Balanced Accuracy', balanced_accuracy),
        ('Accuracy', accuracy),
#         ('F1 Score', f1_score),
        ('MCC', mcc),
        ('Cohen Kappa', cohen_kappa), 
        ('Youdens J', j),
        ('True Positive', TP),
        ('False Negative', FN),
        ('True Negative', TN),
        ('False Positive', FP),
    ])
    ordering_functions = {
        'Sensitivity': lambda x: -x,
        'Specificity': lambda x: -x,
        'Precision': lambda x: -x,
        'Balanced Accuracy': lambda x: -x,
        'Accuracy': lambda x: -x,
#         'F1 Score': lambda x: -x,
        'MCC': lambda x: -x,
        'Cohen Kappa': lambda x: -x,
        'Youdens J': lambda x: -x,
        'True Positive': lambda x: -x,
        'False Negative': lambda x: x,
        'True Negative': lambda x: -x,
        'False Positive': lambda x: x,
    }
    latex_header_conversions = {
        'Sensitivity': 'Sensitivity (TPR)',
        'Specificity': 'Specificity (TNR)',
        'Precision': 'Precision (PPV)',
        'Balanced Accuracy': 'Balanced Accuracy', 
        'Accuracy': 'Accuracy',
#         'F1 Score': 'F1 Score',
        'MCC': 'MCC',
        'Cohen Kappa': 'Cohen Kappa',
        'Youdens J': 'Youdens J',
        'True Positive': 'True Positive',
        'False Negative': 'False Negative',
        'True Negative': 'True Negative',
        'False Positive': 'False Positive',
    }

    # Load submissions data.
    print("Loading submissions...")
    submissions = load_submissions(Stage1Submission, STAGE_1_SUBMISSIONS_DIR_PATH, user_map)
        
    print("Submissions:\n", submissions)
    
    
    # Create submission collection
    print("Generating collection file...")
    
    #define function to generate custom confusion matrix for each submission
    def custom_confusion_matrix(y_true, y_pred, display_labels=None, title=None, output_name=None):
        """
        A function to plot a custom confusion matrix for each submission in Protein-Ligand challenge
        """
        
        #create matrix and flip
        cm = np.flip(confusion_matrix(y_true, y_pred))
        
        #create plot
        fig, ax = plt.subplots(figsize=(8,5))
        axs = sns.heatmap(cm, norm=LogNorm(), annot=True, fmt='d', annot_kws={'fontsize':14})
        axs.plot(ax=ax)
        #cmp = ConfusionMatrixDisplay(cm, display_labels=display_labels)
        #cmp.plot(ax=ax)
        ax.set_title(title, fontsize=14)
        ax.set_xlabel("Predicted", fontsize=14, labelpad=10)
        ax.xaxis.set_ticklabels(['True', 'False'])
        ax.set_ylabel("Experimental", fontsize=14, labelpad=10)
        ax.yaxis.set_ticklabels(['True', 'False'])
        
        return axs


    #Create list of all rankings and sites to iterate over
    rankings = ["Ranked_and_non-ranked", "Ranked", "Non-ranked"]
    sites = ["site-1"]

    #iterates over lists to create statistics for each combination idenpedently
    for ranking in rankings:
        for site in sites:
            OUTPUT_DIRECTORY_PATH_SPECIFIC = '/Users/amezcum1/Desktop/SAMPL9/protein_ligand/Analysis/Analysis-outputs-stage1/{}/{}'.format(ranking,site)
            stage1_submission_collection_specific_file_path = '{}/stage1_submission_collection_{}_{}.csv'.format(OUTPUT_DIRECTORY_PATH_SPECIFIC,ranking,site)

            collection_specific = Stage1SubmissionCollection(submissions,
                                                             experimental_data,
                                                             OUTPUT_DIRECTORY_PATH_SPECIFIC,
                                                             stage1_submission_collection_specific_file_path, 
                                                            ignore_refcalcs=False)

            #Assign False values to non-submitted predictions for molecules. 
            #If a participant only submits "inhibitors" or True values, 
            # everything else should be considered "non-inhibitors" or False. 

            collection_specific.complete_predictions_with_missing_fragments(fragments_file_path=EXPERIMENTAL_DATA_FILE_PATH,
                                                                           submission_collection_file_path=stage1_submission_collection_specific_file_path,
                                                                            ranking=ranking) 

            sns.set_context('talk')
            
            #generate statistics tables for specific collections
            collection_specific.generate_statistics_tables(stats_funcs, subdirectory_path='StatisticsTables',
                                        groupby='SID', site = site,
                                        extra_fields=None, sort_stat='Sensitivity',
                                        ordering_functions=ordering_functions,
                                        latex_header_conversions=latex_header_conversions,
                                        caption='')
            
            #generate confusion matrix plots (flipped)
            #load specific collection into a temp df
            df = pd.read_csv(stage1_submission_collection_specific_file_path)
            
            #make list of submission ids (SIDs)
            sids = []
            for i in df['SID'].unique():
                sids.append(i)
                
            #iterate through submissions and plot confusion matrix
            for sid in sids:
                #make new temp df for each submission
                tmp_df = df.loc[df['SID'] == sid]
                
                #make confusion matrix
                custom_confusion_matrix(tmp_df['Site 1 (exp)'], tmp_df['Site 1 (pred)'],
                                       display_labels=['True','False'], title='SID {} -- {}'.format(sid, tmp_df['Name'].unique().item()))
                
                #save matrix
                plt.savefig('{}/sid-{}_confusion_matrix.png'.format(OUTPUT_DIRECTORY_PATH_SPECIFIC,sid), bbox_inches = 'tight')
            
