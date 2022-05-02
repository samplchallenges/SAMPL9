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
from pkganalysis.stats import (compute_bootstrap_statistics, rmse, mae,
                               me, r2, slope, kendall_tau)


# =============================================================================
# TO DO NOTES
# =============================================================================

# SOME work towards some of these in the submission checking script.
# Specifically, that already correctly loads predictions.

# This script is adapted from SAMPL6/SAMPL7 work by Andrea Rizzi and David Mobley. 

# UPDATES TO PRESENT
# - (DONE) Update host-guest systems, names/etc
# - (DONE) Update weakest and strongest binder data
# - (DONE) Update user-map ( not directly on script)
# - (DONE) Make experimental tables using generate_tables.py (not directly on script)
# - (DONE) Update submission files and give unique method name(s) (not directly on script)
#      - Particularly for participants who included multiple submissions. 
# - (DONE) Update/exclude compounds with previously published values. (WP6 --> WP6-G4 was optional. )
# - (DONE) Update to have "ranked only" and "all" (ranked and non-ranked) analysis output separate
# - (DONE) Include bCD dataset for analysis, re-run analysis.
# - (DONE) Update deprecated ".as_matrix()" function, to ".to_numpy()" 
# - (DONE) Update sids 
# - (DONE) Update stat limits to fit ranked method in plots. 


# TO DO / WIP
# - Update script to generate single plot of correlations of ranked methods for each dataset
#      - fix titles (or change method name), change number of plots, etc. 
# - Update script after figure5 to generate plots of methods that correctly predict strongest/weakest binders. 
# - Update script to commpare similar methods, typically by name (i.e. force field, sampling, charging scheme, water, etc)

# =============================================================================
# CONSTANTS
# =============================================================================

# Paths to input data.
HOST_GUEST_WP6_SUBMISSIONS_DIR_PATH = '/home/amezcum1/SAMPL9/host_guest/Analysis/Submissions/WP6/'
HOST_GUEST_CD_SUBMISSIONS_DIR_PATH = '/home/amezcum1/SAMPL9/host_guest/Analysis/Submissions/CD/'
EXPERIMENTAL_DATA_FILE_PATH = '/home/amezcum1/SAMPL9/experimental_data/experimental_measurements.csv'

# Host color scheme.
HOST_PALETTE = {
    'bCD': '#FFBE0C',
    'HbCD': '#FF0C0C', #First pass
    'CD': '#FFBE0C',
    'WP6': 'C0',
   # 'CD': 'C2',
   # 'bCD': 'C2',
   # 'MGLab_8': 'C2',
   # 'MGLab_9': 'C2',
   # 'MGLab_19': 'C2',
   # 'MGLab_23': 'C2',
   # 'MGLab_24': 'C2',
   # 'MGLab_34': 'C2',
   # 'MGLab_35': 'C2',
   # 'MGLab_36': 'C2',
    'other1': '#1C0F13',  # Licorice
    'other2': '#93A3B1'  # Cadet grey
}


# =============================================================================
# MAIN CHALLENGE HOST-GUEST SUBMISSION
# =============================================================================

class HostGuestSubmission(SamplSubmission):
    """A submission for the main host-guest challenge.

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
    REF_SUBMISSION_SIDS = ['18', '19', '20', '21', '22', '23']  

    # Section of the submission file.
    SECTIONS = {'Predictions', 'Participant name', 'Participant organization', 'Name', 'Software', 'Method', 'Category', 'Ranked'}


    # Sections in CSV format with kwargs to pass to pandas.read_csv().

    CSV_SECTIONS = {
        'Predictions': {'names': ('System ID', '$\Delta$G', 'SEM $\Delta$G', 'd$\Delta$G',
                                  '$\Delta$H', 'SEM $\Delta$H', 'd$\Delta$H'),
                        'index_col': 'System ID'}
    }

    # Acceptable host names (in filenames) and the host IDs they correspond to
    #HOST_NAMES = { 'CLIP': ['CLIP'], 'CD':['bCD', 'MGLab_8', 'MGLab_9', 'MGLab_19', 'MGLab_23', 'MGLab_24', 'MGLab_34', 'MGLab_35', 'MGLab_36'],
    #            'GDCC':['OA', 'exoOA'] }

    HOST_NAMES = { 'WP6': ['WP6'], 'CD': ['bCD', 'HbCD']}
    #HOST_NAMES = { 'WP6': ['WP6']}

    RENAME_METHODS = {}

    def __init__(self, file_path, user_map):
        super().__init__(file_path, user_map)

        file_name = os.path.splitext(os.path.basename(file_path))[0]
        self.file_name = file_name

        #TO DO:  Not sure if I'm going to use the immediately following for anything
        file_name_simple = file_name.replace('_','-')
        file_data = file_name_simple.split('-')
        self.host_name = file_data[0]

        # Load predictions.
        sections = self._load_sections(file_path)  # From parent-class.
        self.data = sections['Predictions']  # This is a list
        self.data = pd.DataFrame(data=self.data) # Now a DataFrame
        try:
            self.name = self.RENAME_METHODS[sections['Name'][0]]
        except KeyError:
            self.name = sections['Name'][0]


        # Add host name column to predictions.
        self.host_name = file_data[0].upper()
        self.data['host_name'] = self.host_name
        assert self.host_name in self.HOST_NAMES

        # Store participant name, organization, method category
        self.participant = sections['Participant name'][0].strip()
        self.category = sections['Category'][0].strip()
        self.organization = sections['Participant organization'][0].strip()
        self.ranked = sections['Ranked'][0].strip() =='True'

        # Required system System IDs. For WP6 G4 is not required. 
        WP6_guests = ['G1', 'G2', 'G3', 'G5', 'G6', 'G7', 'G8', 'G9', 'G10', 'G11', 'G12', 'G13']
        CD_guests = ['TDZ','TFP', 'PMZ', 'PMT', 'CPZ']

        self.REQUIRED_SYSTEM_IDs = {'WP6':[f'WP6-{guest}' for guest in WP6_guests], 
                                    'CD':['bCD-TDZ', 'bCD-TFP', 'bCD-PMZ', 'bCD-PMT', 'b-CD-CPZ'] + ['HbCD-TDZ', 'HbCD-TFP','HbCD-PMZ', 'HbCD-PMT', 'HbCD-CPZ']}
        #self.REQUIRED_SYSTEM_IDs = {'WP6':[f'WP6-{guest}' for guest in WP6_guests]}

        # Check if this is a test submission.
        if self.sid in self.TEST_SUBMISSION_SIDS:
            raise IgnoredSubmissionError('This submission has been used for tests.')

        # Check if this is a reference submission
        self.reference_submission = False
        if self.sid in self.REF_SUBMISSION_SIDS:
            self.reference_submission = True

    def __add__(self, other):
        """Merge the data of the two submission."""
        merged_submission = copy.deepcopy(self)
        merged_submission.sid = '{} + {}'.format(*sorted([self.sid, other.sid]))
        merged_submission.host_name = '{} + {}'.format(*sorted([self.host_name, other.host_name]))

        # Check if this is already a merged submission.
        if isinstance(merged_submission.file_name, list):
            merged_submission.file_name = sorted([*merged_submission.file_name, other.file_name])
        else:
            merged_submission.file_name = sorted([merged_submission.file_name, other.file_name])
        merged_submission.data = pd.concat([merged_submission.data, other.data])
        return merged_submission

    def split(self, names_to_separate):
        """Take a host-guest submission that spans multiple hosts (with system IDs including host name and guest name), and split it into multiple submissions
        which have the same metadata but only contain the data for the individual hosts. The resulting submissions have updated `host_name` fields also.

        Takes a list of host names (as used for the individual data points) to separate based on.

        Returns a list of the new HostGuestSubmission objects, of length equal to `names_to_separate`"""

        # Find how many submissions we're making and make new submissions, duplicating old
        n_submissions = len(names_to_separate)
        new_submissions = [copy.deepcopy(self) for i in range(n_submissions)]



        # Build list of system IDs we want
        for (idx, submission) in enumerate(new_submissions):
            desired_IDs = []
            for system_id, series in submission.data[['$\Delta$G', 'd$\Delta$G', '$\Delta$H']].iterrows():
                tmp = system_id.split('-')
                if tmp[0] == names_to_separate[idx]:
                    desired_IDs.append( system_id )

            # Grab just that data and store
            new_submissions[idx].data = submission.data.loc[desired_IDs]
            # Change the host name to what's correct for this host
            new_submissions[idx].data.host_name = names_to_separate[idx]
            new_submissions[idx].host_name = names_to_separate[idx]

        return new_submissions



# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

class HostGuestSubmissionCollection:
    """A collection of HostGuestSubmissions."""

    FREE_ENERGY_CORRELATION_PLOT_DIR = 'FreeEnergyCorrelationPlots'
    ENTHALPIES_CORRELATION_PLOT_DIR = 'EnthalpiesCorrelationPlots'
    MOLECULE_CORRELATION_PLOT_PATH = 'molecules_error.pdf'

    _ROW_HEIGHT = 0.25


    def __init__(self, submissions, experimental_data, output_directory_path, ignore_refcalcs = True, ranked_only = True, allow_multiple = False):
        # Use "allow_multiple" if we're using a submission collection which is aggregated across several hosts
        # in which case participants may have a ranked submission in each separate host.

        # Build full free energy table.
        data = []

        # Participant names we've found so far; tracked to ensure no one has more than one
        # ranked submission
        self.participant_names_ranked = []

        # Submissions free energies and enthalpies.
        for submission in submissions:
            # Ignore reference calculations, if applicable
            if submission.reference_submission and ignore_refcalcs:
                continue

            if ranked_only and not submission.ranked:
                continue

            # Store names associated with ranked submission, skip if they submitted multiple (only if we need to check for duplicate authors)
            if submission.ranked and not allow_multiple:
                if not submission.participant in self.participant_names_ranked:
                    self.participant_names_ranked.append(submission.participant)
                else:
                    print(f"Error: {submission.participant} submitted multiple ranked submissions.")
                    continue

            for system_id, series in submission.data[['$\Delta$G', 'd$\Delta$G', '$\Delta$H']].iterrows():
                try: #Load experimental data
                    free_energy_expt = experimental_data.loc[system_id, '$\Delta$G']
                    enthalpy_expt = experimental_data.loc[system_id, '$\Delta$H']
                except KeyError: #But if there's no data for this one, skip
                    continue
                free_energy_calc = series['$\Delta$G']
                free_energy_calc_sem = series['d$\Delta$G']
                enthalpy_calc = series['$\Delta$H']
                ranked = submission.ranked


                data.append({
                    'sid': submission.sid,
                    'participant': submission.participant,
                    'name': submission.name,
                    #'method': self._assign_paper_method_name(submission.name),
                    'method': submission.name, # Make this duplicate name for now, as right now name does somewhat describe method. TO DO
                    'system_id': system_id,
                    'host_name': submission.host_name,
                    '$\Delta$G (calc) [kcal/mol]': free_energy_calc,
                    'd$\Delta$G (calc) [kcal/mol]': free_energy_calc_sem,
                    '$\Delta$G (expt) [kcal/mol]': free_energy_expt,
                    '$\Delta\Delta$G error (calc - expt)  [kcal/mol]': free_energy_calc - free_energy_expt,
                    '$\Delta$H (calc) [kcal/mol]': enthalpy_calc,
                    '$\Delta$H (expt) [kcal/mol]': enthalpy_expt,
                    '$\Delta\Delta$H error (calc - expt)  [kcal/mol]': enthalpy_calc - enthalpy_expt
                })

        # Transform into Pandas DataFrame.
        self.data = pd.DataFrame(data=data)
        self.output_directory_path = output_directory_path

        # Create general output directory.
        os.makedirs(self.output_directory_path, exist_ok=True)

    @staticmethod
    def _assign_method_class(name):
        # Roughly groups similar methods.
        method_classes = {
            'AMOEBA/BAR/Tinker': 'Alchemical/AMOEBA',
            'FS-DAM/GAFF2/TIP3P': 'Alchemical/Nonequilibrium'
        }
        if name in method_classes:
            return method_classes[name]
        if name.startswith('BSSE') or name.startswith('DFT') or name.startswith('SQM'):
            return 'QM'
        if 'MMPBSA' in name:
            return 'MMPBSA'
        if 'Umbrella Sampling' in name or (name.startswith('US') and name[-1] == '1'):
            return 'Umbrella Sampling'
        if name.startswith('US'):
            n_submission = int(name.split('_')[1])
            if n_submission == 1:
                return 'Umbrella Sampling'
            elif n_submission == 2:
                return 'Umbrella Sampling/Fitted'
            elif n_submission in {3, 5, 8, 10, 13, 15, 18, 20}:
                return 'Movable Type'
            else:
                return 'Movable Type/Fitted'
        if 'Force-Matching' in name:
            return 'Alchemical/Force-Matching'
        if name.startswith('SOMD'):
            if name[-1] == 'D':
                return 'Alchemical/GAFF/Fitted'
            else:
                return 'Alchemical/GAFF'
        if 'GAFF' in name:
            return 'Alchemical/GAFF'
        if name.startswith('FEP'):
            return 'Alchemical/Relative'
        return name

    # TO DO: The following function does not pertain to this challenge/needs updating if we even retain
    @staticmethod
    def _assign_paper_method_name(name):
        return name
        # Convert from submission method name to the name used in the paper.
    #    method_names = {
    #        'DDM/GAFF/AM1-BCC/TIP3P': 'DDM-GAFF',
    #        'HREM/BAR/RESP/Force-Matching/TIP3P': 'DDM-FM',
    #        'DDM/Force-Matching/FEP/HREM/MBAR': 'DDM-FM-QMMM',
    #        'BSSE-corrected RI-B3PW91 (SMD)/CBS': 'DFT(B3PW91)',
    #        'BSSE-corrected RI-B3PW91-D3 (SMD)/CBS': 'DFT(B3PW91)-D3',
    #        'DFT-opt': 'DFT(TPSS)-D3',
    #        'FEP-MM': 'RFEC-GAFF2',
    #        'FEP-QM/MM': 'RFEC-QMMM',
    #        'FS-DAM/GAFF2/TIP3P': 'FSDAM',
    #        'EKEN-DIAZ/MD/MMPBSA': 'MMPBSA-GAFF',
    #        'SQM-opt': 'SQM(PM6-DH+)',
    #        'AMOEBA/BAR/Tinker': 'DDM-AMOEBA',
    #        'Umbrella Sampling/TIP3P': 'US-CGenFF',
    #        'US/PMF/MT/MD_1': 'US-GAFF',
    #        'US/PMF/MT/MD_2': 'US-GAFF-C',
    #    }
    #    try:
    #        return method_names[name]
    #    except KeyError:
    #        pass
#
#        if name.startswith('SOMD/AM1BCC-GAFF-TIP3P'):
#            paper_name = 'SOMD-' + name[-1]
#            if 'NOBUFFER' in name:
#                paper_name += '-nobuffer'
#        if name.startswith('US/PMF/MT/MD'):
#            submission_number = int(name.rsplit('_', 1)[-1])
#            identifier = ''
#            # Potential
#            if 18 <= submission_number <= 19:
#                identifier += 'K'
#            else:
#             # Input structure identifier
#            if 3 <= submission_number <= 7:  # MD_relaxed
#                identifier += 'D'
#            elif 8 <= submission_number <= 12:  # minimum from Umbrella Sampling
#                identifier += 'U'
#            elif 13 <= submission_number <= 19 or 26 <= submission_number <= 27:  # MT_minimum
#                identifier += 'T'
#            # States
#            if submission_number in [3, 4, 8, 9, 13, 14, 20, 21, 22, 23, 24, 25, 26, 27]:
#                identifier += '3'
#            else:
#                identifier += '1'
#            # Correction
#            if submission_number in [3, 5, 8, 10, 13, 15, 18, 20]:  # No
#                identifier += 'N'
#            elif submission_number in [4, 6, 9, 11, 14, 16, 19, 21]:  # Yes (linear)
#                identifier += 'L'
#            elif submission_number in [7, 12, 17, 22]:  # Yes only intercept
#                identifier += 'O'
#            elif submission_number in [23, 26]:  # Yes^1
#                identifier += 'U'
#            elif submission_number in [27, 24]:  # Yes^2
#                identifier += 'S'
#            else:  # Yes^2 only intercept
#                identifier += 'Z'
#            paper_name = 'MovTyp-' + identifier
#        if 'NULL' in name:
#            paper_name = name
#        return paper_name

    def generate_correlation_plots(self):
        # Free energy correlation plots.
        self._generate_correlation_plots(x='$\Delta$G (expt) [kcal/mol]', y='$\Delta$G (calc) [kcal/mol]',
                                         directory_path=self.FREE_ENERGY_CORRELATION_PLOT_DIR)
        # Enthalpies correlation plots.
        self._generate_correlation_plots(x='$\Delta$H (expt) [kcal/mol]', y='$\Delta$H (calc) [kcal/mol]',
                                         directory_path=self.ENTHALPIES_CORRELATION_PLOT_DIR)

    def _generate_correlation_plots(self, x, y, directory_path):
        output_dir_path = os.path.join(self.output_directory_path, directory_path)
        os.makedirs(output_dir_path, exist_ok=True)
        for sid in self.data.sid.unique():
            data = self.data[self.data.sid == sid]

            # If this is a merged submission, we need a hue.
            host_names = data.host_name.unique()
            if len(host_names) > 1:
                hue = 'host_name'
                title = '{} ({})'.format(data.name.unique()[0], sid)
            else:
                hue = None
                title = '{} - {} ({})'.format(data.name.unique()[0], host_names[0], sid)

            # Check if enthalpies were computed.
            if data[y].isnull().any():
                continue

            plt.close('all')
            plot_correlation(x=x, y=y, data=data, title=title, hue=hue)
            plt.tight_layout(pad=0.2)
            # plt.show()
            output_path = os.path.join(output_dir_path, '{}.pdf'.format(sid))
            plt.savefig(output_path)


    def generate_molecules_plot(self):
        # Correlation plot by molecules.
        plt.close('all')
        n_rows = len(self.data.system_id.unique())
        fig, ax = plt.subplots(figsize=(6, 0.4*n_rows))
        sns.violinplot(y='system_id', x='$\Delta\Delta$G error (calc - expt)  [kcal/mol]',
                       data=self.data, linewidth=1.0, inner='point', cut=0, ax=ax)
        plt.tight_layout(pad=0.2)
        # plt.show()
        plt.savefig(os.path.join(self.output_directory_path, self.MOLECULE_CORRELATION_PLOT_PATH))

    def generate_statistics_tables(self, stats_funcs, subdirectory_path, groupby,
                                   extra_fields=None, sort_stat=None,
                                   ordering_functions=None, latex_header_conversions=None,
                                   caption='', ignore_refcalcs = True):
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
            return s.replace('_', '\_')

        extra_fields_latex = [escape(extra_field) for extra_field in extra_fields]

        file_base_name = 'statistics'
        directory_path = os.path.join(self.output_directory_path, subdirectory_path)

        stats_names, stats_funcs = zip(*stats_funcs.items())
        ci_suffixes = ('', '_lower_bound', '_upper_bound')

        # Compute or read the bootstrap statistics from the cache.
        cache_file_path = os.path.join(self.output_directory_path, 'bootstrap_distributions.p')
        all_bootstrap_statistics = self._get_bootstrap_statistics(groupby, stats_names, stats_funcs,
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
        print()

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

        # Reorder columns that were scrambled by going through a dictionaries.
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

    def plot_bootstrap_distributions(self, stats_funcs, subdirectory_path, groupby,
                                     ordering_functions=None, latex_header_conversions=None,
                                     stats_limits=None, exclusions=frozenset(),
                                     shaded=frozenset(), figure_width=7.25, output_file_suffix='',
                                     **violinplot_kwargs):
        """Generate a violin plot for the bootstrap distribution of the statistics.

        Parameters
        ----------
        subdirectory_path : str
            The path of the output directory relative to self.output_directory_path.
        groupby : str
            The name of the data column to be used to compute the statistics.
            For example, 'name' to obtain statistics about individual methods,
            'system_id' to compute statistics by molecules.
        ordering_functions : dict
            Dictionary statistic_name -> ordering_function(stats), where
            ordering_function determines how to rank the the groups by
            statistics.
        exclusions : set
            Which groups to exclude from the plot.
        figure_width : float, optional
            The width of the figure for each statistic (default is 7.25).
        output_file_suffix : str, optional
            A suffix including the file extension. The image will be saved with
            the name statisticname_bootstrap_distributions + output_file_suffix.
            By default, no suffix is used and the image is saved in PDF format.

        """
        if stats_limits is None:
            stats_limits = {}

        directory_path = os.path.join(self.output_directory_path, subdirectory_path)
        stats_names, stats_funcs = zip(*stats_funcs.items())

        # Compute or read the bootstrap statistics from the cache.
        ordering_data, statistics_plot = self._get_bootstrap_distribution_plot_data(
            groupby, stats_names, stats_funcs)

        # Create CSV and JSON tables (correct LaTex syntax in column names).
        os.makedirs(directory_path, exist_ok=True)

        # Violin plots by statistics across submissions.
        n_groups = len(ordering_data.index) - len(exclusions)
        for stats_name in stats_names:
            plt.close('all')
            fig, ax = plt.subplots(figsize=(figure_width, self._ROW_HEIGHT*(n_groups + 1)))

            stats_name_latex = latex_header_conversions.get(stats_name, stats_name)
            data = statistics_plot[statistics_plot.statistics_name == stats_name_latex]

            # Determine the order in which to display groups.
            ordering_function = ordering_functions.get(stats_name, lambda x: x)
            order = sorted(ordering_data[stats_name].items(), key=lambda x: ordering_function(x[1]))
            order = [group for group, stats in order if group not in exclusions]

            # Plot boot strap distributions.
            sns.violinplot(x='value', y='ID', data=data, linewidth=0.8,
                           order=order, scale='width', width=2.8*self._ROW_HEIGHT,
                           cut=0, ax=ax, **violinplot_kwargs)

            self._modify_violinplot(ax, stats_name)

            # Configure axes.
            if stats_name in stats_limits:
                ax.set_xlim(stats_limits[stats_name])
            ax.set_xlabel(stats_name_latex)
            ax.set_ylabel('')

            # Create shaded area.
            for order_idx, group in enumerate(order):
                if group in shaded:
                    x_limits = ax.get_xlim()
                    y = np.full(len(x_limits), float(order_idx))
                    # Set the zorder of the violinplot element so that the shaded area remains below.
                    ax.fill_between(x_limits, y - 0.5, y + 0.5, alpha=0.1, color='0.1', linewidth=0.0, zorder=0.5)
                    # Make sure the shaded area doesn't expand the old axis limits.
                    ax.set_xlim(x_limits)

            # Remove legend if present.
            if ax.legend_ is not None:
                ax.legend_.remove()
            # Move ticks closer to axes.
            ax.tick_params(pad=2.0)
            plt.tight_layout(pad=0.24)
            # plt.tight_layout()
            # Make sure that here are enough ticks.
            if len(ax.get_xticks()) < 5 and stats_name in stats_limits:
                stats_lim = stats_limits[stats_name]
                ax.xaxis.set_ticks(np.linspace(stats_lim[0], stats_lim[1], num=6))

            # plt.show()
            # Check if we need to use a default image format.
            if output_file_suffix == '':
                suffix = '.pdf'
            else:
                suffix = '_' + output_file_suffix
            output_file_path = os.path.join(directory_path, stats_name) + '_bootstrap_distributions' + suffix
            plt.savefig(output_file_path)

    def _get_bootstrap_distribution_plot_data(self, groupby, stats_names, stats_funcs):
        """Return the dataframes with the statistics necessary to plot.

        Returns
        -------
        ordering_data : pandas.Dataframe
            A dataframe containing all the statistics that can be used to decide
            the order of the bootstrap distributions to plot.
        statistics_plot_data : pandas.Dataframe
            A dataframe containing all the bootstrap samples.
        """
        cache_file_path = os.path.join(self.output_directory_path, 'bootstrap_distributions.p')
        all_bootstrap_statistics = self._get_bootstrap_statistics(groupby, stats_names, stats_funcs,
                                                                  cache_file_path=cache_file_path)

        # Collect the records for the DataFrames.
        ordering_data = []
        statistics_plot = []

        groups = self.data[groupby].unique()
        for i, group in enumerate(groups):
            print('\rCollecting bootstrap statistics for {} {} ({}/{})'
                  ''.format(groupby, group, i+1, len(groups)), end='')

            # Isolate bootstrap statistics.
            bootstrap_statistics = all_bootstrap_statistics[group]

            ordering_data_record = {}
            for stats_name, (stats, (lower_bound, upper_bound), bootstrap_samples) in bootstrap_statistics.items():
                stats_name_latex = latex_header_conversions.get(stats_name, stats_name)
                # Keep the data used for the order.
                # TODO check primary_hue
                ordering_data_record[stats_name] = stats
                # For the violin plot, we need all the bootstrap statistics series.
                for bootstrap_sample in bootstrap_samples:
                    statistics_plot.append(dict(ID=group, statistics_name=stats_name_latex,
                                                value=bootstrap_sample))
            ordering_data.append(dict(ID=group, **ordering_data_record))
        print()

        ordering_data = pd.DataFrame(ordering_data)
        ordering_data.set_index('ID', inplace=True)
        statistics_plot = pd.DataFrame(statistics_plot)
        return ordering_data, statistics_plot

    def _get_bootstrap_statistics(self, groupby, stats_names, stats_funcs, cache_file_path):
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
                            all_bootstrap_statistics[group][stats_name] = cached_bootstrap_statistics[method_name][stats_name]
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
                  ''.format(groupby, group, i+1, len(groups)), end='')

            # Select the group data.
            data = self.data[self.data[groupby] == group]

            # Check if SEMs for the free energies are reported.
            sems = data['d$\Delta$G (calc) [kcal/mol]'].values
            if np.any(np.isnan(sems)):
                sems = None
            else:  # Add a column of SEMs = 0.0 for the experimental values.
                sems = np.array([(0.0, sem) for sem in sems])

            # Compute bootstrap statistics.
            data = data[['$\Delta$G (expt) [kcal/mol]', '$\Delta$G (calc) [kcal/mol]']]
            new_bootstrap_statistics = compute_bootstrap_statistics(data.to_numpy(), group_stats_funcs, sems=sems,
                                                                    n_bootstrap_samples=10000)

            # Update the returned value with the statistics just computed.
            new_boostrap_statistics = {group_stats_names[i]: new_bootstrap_statistics[i]
                                       for i in range(len(group_stats_funcs))}
            group_bootstrap_statistics.update(new_boostrap_statistics)

        # Cache the computed statistics on disk. Create output directory if necessary.
        if cache_updated:
            os.makedirs(os.path.dirname(cache_file_path), exist_ok=True)
            with open(cache_file_path, 'wb') as f:
                pickle.dump(all_bootstrap_statistics, f)

        return all_bootstrap_statistics

    def _modify_violinplot(self, ax, stats_name):
        pass


# =============================================================================
# MERGE SUBMISSIONS AND COLLECTIONS
# =============================================================================

def merge_submissions(submissions, discard_not_matched=True):
    # Find all host names.
    host_names = set([submission.host_name for submission in submissions])

    # Find submissions that have the same name.
    submissions_by_name = {}
    for submission in submissions:
        try:
            submissions_by_name[submission.name].append(submission)
        except KeyError:
            submissions_by_name[submission.name] = [submission]

    # Merge TEMOA/TEETOA submissions that use the same method into a single submission object.
    merged_submissions = []
    for method_name, method_submissions in submissions_by_name.items():
        # Check that the submissions come from the same participant,
        # and that it's the same method applied to different systems.
        assert len(method_submissions) <= len(host_names)
        assert len(set([submission.participant for submission in method_submissions])) == 1
        assert len(set([submission.host_name for submission in method_submissions])) == len(method_submissions)

        # Discard methods that were run on only a subset of the hosts.
        if len(method_submissions) == len(host_names) or not discard_not_matched:
            if len(method_submissions) == 1:
                merged_submissions.append(method_submissions[0])
            else:
                merged_submissions.append(sum(method_submissions[1:], method_submissions[0]))

    return merged_submissions


# =============================================================================
# UTILITIES TO GENERATE FIGURES FOR THE PAPER
# =============================================================================

class SplitBootstrapSubmissionCollection(HostGuestSubmissionCollection):
    """Two collections merged, which allows to plot split bootstrap distributions.

    The violin plots are ordered according to collection1.

    Parameters
    ----------
    collection1 : HostGuestSubmissionCollection
        The first collection.
    collection2 : HostGuestSubmissionCollection
        The second collection.
    hue : str
        The name of the hue parameter in seaborn violinplot.
    collection1_hue : str
        The hue value for collection1.
    collection2_hue : str
        The hue value for collection2.
    output_directory_path : str
        The path of the directory where to save the results of the analysis.
    """

    def __init__(self, collection1, collection2, hue, collection1_hue, collection2_hue, output_directory_path):
        self.hue = hue
        self.collections = collections.OrderedDict([
            (collection1_hue, collection1),
            (collection2_hue, collection2)
        ])

        # Create output directory if needed.
        self.output_directory_path = output_directory_path
        os.makedirs(self.output_directory_path, exist_ok=True)

    def generate_paper_table(self, stats_funcs, exclusions):
        # Call _get_bootstrap_distribution_plot_data() to create self._collection_statistics.
        stats_names, stats_funcs = zip(*stats_funcs.items())
        self._get_bootstrap_distribution_plot_data('method', stats_names, stats_funcs)
        data = self._collections_statistics  # Short-cut.

        # Remove exclusions.
        data = data[~data.ID.isin(exclusions)]

        # We'll print all the methods in alphabetical order.
        methods = sorted(data.ID.unique())

        # Build a table with statistics of TEMOA/TEETOA and CB8 side by side.
        # The first column of the table list the methods.
        table = [['{:15s}'.format(method)] for method in methods]
        for hue in self.collections:
            for method_idx, method in enumerate(methods):
                for stats_name in stats_names:
                    try:
                        stats = data[(data[self.hue] == hue) & (data.ID == method)]
                        stats = '{:.1f} ({:.1f}) & [{:.1f}, {:.1f}]'.format(
                            stats[stats_name].values[0], stats[stats_name + '_mean'].values[0],
                            stats[stats_name + '_low'].values[0], stats[stats_name + '_up'].values[0])
                    except IndexError:
                        # No data for this combination of dataset and method.
                        stats = ' & '
                    table[method_idx].append(stats)

        # Plot table.
        table_str = ''
        for row in table:
            table_str += ' & '.join(row) + ' \\\\\n'
        with open('temp.txt', 'w') as f:
            f.write(table_str)

    def _get_bootstrap_distribution_plot_data(self, groupby, stats_names, stats_funcs):
        """Return the dataframes with the statistics necessary to plot.

        Returns
        -------
        ordering_data : pandas.Dataframe
            A dataframe containing all the statistics that can be used to decide
            the order of the bootstrap distributions to plot.
        statistics_plot_data : pandas.Dataframe
            A dataframe containing all the bootstrap samples.
        """
        # Collect the records for the DataFrames.
        statistics_plot = []
        # Cache the statistics of the collections as we'll need them in _modify_violinplot().
        self._collections_statistics = []

        # Gather data for both collection distributions.
        for collection_hue, collection in self.collections.items():
            cache_file_path = os.path.join(collection.output_directory_path, 'bootstrap_distributions.p')
            all_bootstrap_statistics = collection._get_bootstrap_statistics(groupby, stats_names, stats_funcs,
                                                                            cache_file_path=cache_file_path)

            groups = collection.data[groupby].unique()

            for group_idx, group in enumerate(groups):
                print('\rCollecting bootstrap statistics for {} {} ({}/{})'
                      ''.format(groupby, group, group_idx+1, len(groups)), end='')
                # Isolate bootstrap statistics.
                bootstrap_statistics = all_bootstrap_statistics[group]

                ordering_data_record = {}
                for stats_name, (stats, (lower_bound, upper_bound), bootstrap_samples) in bootstrap_statistics.items():
                    stats_name_latex = latex_header_conversions.get(stats_name, stats_name)
                    # Keep the data used for the order.
                    ordering_data_record[stats_name] = stats
                    # These are used to create the table for the paper.
                    ordering_data_record[stats_name + '_low'] = lower_bound
                    ordering_data_record[stats_name + '_up'] = upper_bound
                    ordering_data_record[stats_name + '_mean'] = bootstrap_samples.mean()
                    # For the violin plot, we need all the bootstrap statistics series.
                    for bootstrap_sample in bootstrap_samples:
                        statistics_plot.append(dict(ID=group, statistics_name=stats_name_latex,
                                                    value=bootstrap_sample, **{self.hue: collection_hue}))
                # Save actual statitics.
                self._collections_statistics.append(dict(ID=group, **ordering_data_record,
                                                         **{self.hue: collection_hue}))

            del all_bootstrap_statistics  # This can be huge when the number of bootstrap cycles is high.
            print()

        statistics_plot = pd.DataFrame(statistics_plot)
        self._collections_statistics = pd.DataFrame(self._collections_statistics)

        # # Print mean statistics.
        # temp = self._collections_statistics
        # for hue in ['TEMOA/TEETOA', 'CB8']:
        #     for stats_name in stats_names:
        #         print(hue, stats_name, temp[temp[self.hue] == hue][stats_name].mean())

        # We order using collection1 statistics (unless the group has only collection2).
        (collection1_hue, collection1), (collection2_hue, collection2) = self.collections.items()
        ordering_data = self._collections_statistics[self._collections_statistics[self.hue] == collection1_hue]
        groups2_unique = set(collection2.data[groupby].unique()) - set(collection1.data[groupby].unique())
        if len(groups2_unique) != 0:
            groups2_data = self._collections_statistics[(self._collections_statistics[self.hue] == collection2_hue) &
                                                        (self._collections_statistics.ID.isin(groups2_unique))]
            ordering_data = pd.concat([ordering_data, groups2_data], ignore_index=True)
        ordering_data.set_index('ID', inplace=True)

        return ordering_data, statistics_plot

    def _modify_violinplot(self, ax, stats_name):
        """Add the real statistics to the violin plot."""
        for group_idx, group in enumerate(ax.get_yticklabels()):
            for collection_idx, (collection_hue, collection) in enumerate(self.collections.items()):
                collection_data = self._collections_statistics[(self._collections_statistics.ID == group.get_text()) &
                                                               (self._collections_statistics[self.hue] == collection_hue)]
                try:
                    collection_stats = collection_data[stats_name].values[0]
                except IndexError:
                    # The group doesn't have a statistics for the hue.
                    continue
                sign = -1 if collection_idx == 0 else 1
                y = np.linspace(group_idx + sign * self._ROW_HEIGHT*0.5, group_idx + sign * self._ROW_HEIGHT * 1.0)
                x = [collection_stats for _ in y]
                ax.plot(x, y, c='white', lw=2.0,)
                ax.plot(x, y, c='black', lw=1.8, alpha=0.85)


# =============================================================================
# MAIN
# =============================================================================

if __name__ == '__main__':
    # Read experimental data.
    with open(EXPERIMENTAL_DATA_FILE_PATH, 'r') as f:
        # experimental_data = pd.read_json(f, orient='index')
        #names = ('System ID', 'name', 'SMILES', '$Kd_1$', 'd$Kd_1$','$\Delta H_1$', 'd$\Delta H_1$', '$Kd_2$', 'd$Kd_2$','$\Delta H_2$', 'd$\Delta H_2$',
        #         'T$\Delta$S', 'dT$\Delta$S', 'n', '$Ka_1$', 'd$Ka_1$', '$Ka_2$', 'd$Ka_2$', '$Ka$','d$Ka','$\Delta$H', 'd$\Delta$H',
        #         '$\Delta$G', 'd$\Delta$G')
        #names = ('System ID', 'name', 'SMILES', '$Ka_1$', 'd$Ka_1$', '$Ka_2$', 'd$Ka_2$','$\Delta H_1$', 'd$\Delta H1$', 
        #        '$\Delta H_2$', 'd$\Delta H_2$', 'T$\Delta$S', 'dT$\Delta$S', 'n', '$Ka$','d$Ka','$\Delta$H', 
        #        'd$\Delta$H', '$\Delta$G', 'd$\Delta$G')
        names = ('System ID', 'name', 'SMILES', '$\Delta$G', 'd$\Delta$G', '$\Delta$H', 'd$\Delta$H', 'T$\Delta$S', 
                'dT$\Delta$S', 'n', '$Ka$', 'd$Ka$')

        experimental_data = pd.read_csv(f, sep=';', names=names, index_col='System ID', skiprows=1)
        #Don't read experimental values which are non-numeric -- e.g. the experiments didn't work/weren't done
        #experimental_data = experimental_data.dropna(subset=['$\Delta$G', '$Ka$', 'd$Ka$', '$\Delta$H',
        #    'd$\Delta$H', 'd$\Delta$G', 'T$\Delta$S', 'dT$\Delta$S'])
        experimental_data = experimental_data.dropna(subset=['$\Delta$G'])

    # Convert numeric values to dtype float.
    for col in experimental_data.columns[3:]:
        experimental_data[col] = pd.to_numeric(experimental_data[col], errors='coerce')

    # Rename CB8-G12a to CB8-G12 since we'll print only this. PERTAINS TO SAMPL6
    #id_index = np.where(experimental_data.index.values == 'CB8-G12a')[0][0]
    #experimental_data.index.values[id_index] = 'CB8-G12'

    # Import user map.
    try:
        with open('/home/amezcum1/SAMPL9/host_guest/Analysis/SAMPL9-user-map-HG.csv', 'r') as f:
            user_map = pd.read_csv(f)
    except FileNotFoundError:
        user_map=None
        print("Warning: No user map found.")

    # Configuration: statistics to compute.
    stats_funcs = collections.OrderedDict([
        ('RMSE', rmse),
        ('MAE', mae),
        ('ME', me),
        ('R2', r2),
        ('m', slope),
        ('kendall_tau', kendall_tau)
    ])
    ordering_functions = {
        'ME': lambda x: abs(x),
        'R2': lambda x: -x,
        'm': lambda x: abs(1 - x),
        'kendall_tau': lambda x: -x
    }
    latex_header_conversions = {
        'R2': 'R$^2$',
        'RMSE': 'RMSE [kcal/mol]',
        'MAE': 'MAE [kcal/mol]',
        'ME': 'ME [kcal/mol]',
        'kendall_tau': '$\\tau$',
    }
    stats_limits = {
        'RMSE': (0, 25.0),
        'MAE': (0, 25),
        'ME': (-20, 20),
        'R2': (0, 1),
        'm': (-5, 5),
        'kendall_tau': (-1, 1),
        } 

    # Statistics by molecule.
    stats_funcs_molecules = collections.OrderedDict([
        ('RMSE', rmse),
        ('MAE', mae),
        ('ME', me),
    ])

    # Instantiate two sets of HostGuestSubmissionCollections -- one which ignores reference calculations
    # and one which doesn't.

    # Load submissions data. For now only CB8 and GDCC
    print("Loading WP6 submissions")
    submissions_wp6 = load_submissions(HostGuestSubmission, HOST_GUEST_WP6_SUBMISSIONS_DIR_PATH, user_map)
    print("Loading CD submissions")
    submissions_cd = load_submissions(HostGuestSubmission, HOST_GUEST_CD_SUBMISSIONS_DIR_PATH, user_map)

    # Make split submissions for the two hosts within GDCC
    submissions_bcd = []
    submissions_hbcd = []
    for submission in submissions_cd:
        a, b = submission.split(['bCD', 'HbCD'])
        submissions_bcd.append(a)
        submissions_hbcd.append(b)

    # Make a set of all the submissions
    submissions_all = submissions_wp6 + submissions_bcd + submissions_hbcd
    #submissions_all = submissions_bcd + submissions_hbcd

    # Make directories for output
    if not os.path.isdir('../Ranked_Accuracy'): os.mkdir('../Ranked_Accuracy')
    if not os.path.isdir('../Ranked_Accuracy/MoleculesStatistics'): os.mkdir('../Ranked_Accuracy/MoleculesStatistics')
    if not os.path.isdir('../Ranked_Accuracy/PaperImages'): os.mkdir('../Ranked_Accuracy/PaperImages')

    if not os.path.isdir('../All_Accuracy'): os.mkdir('../All_Accuracy')
    if not os.path.isdir('../All_Accuracy/MoleculesStatistics'): os.mkdir('../All_Accuracy/MoleculesStatistics')
    if not os.path.isdir('../All_Accuracy/PaperImages'): os.mkdir('../All_Accuracy/PaperImages')

    # Create submission collections
    print("Creating submission collection for WP6")
    collection_wp6 = HostGuestSubmissionCollection(submissions_wp6, experimental_data,
                                                  output_directory_path='../Ranked_Accuracy/WP6')

    print("Creating submission collection for all WP6, including non-ranked")
    collection_wp6_nonranked = HostGuestSubmissionCollection(submissions_wp6, experimental_data,
            output_directory_path='../All_Accuracy/WP6', 
            ignore_refcalcs = False, ranked_only = False)

    print("Creating submission collection for CD collectively")
    collection_cd = HostGuestSubmissionCollection(submissions_cd, experimental_data,
                                                   output_directory_path='../Ranked_Accuracy/CD')

    print("Creating submission collection for all CD collectively, including non-ranked")
    collection_cd_nonranked = HostGuestSubmissionCollection(submissions_cd, experimental_data,
            output_directory_path='../All_Accuracy/CD',
            ignore_refcalcs = False, ranked_only = False)

    print("Creating submission collection for bCD")
    collection_bcd = HostGuestSubmissionCollection(submissions_bcd, experimental_data,
                                                    output_directory_path='../Ranked_Accuracy/bCD')

    print("Creating submission collection for bCD, including non-ranked")
    collection_bcd_nonranked = HostGuestSubmissionCollection(submissions_bcd, experimental_data,
            output_directory_path='../All_Accuracy/bCD',
            ignore_refcalcs = False, ranked_only = False)

    print("Creating submission collection for HbCD")
    collection_hbcd = HostGuestSubmissionCollection(submissions_hbcd, experimental_data,
                                                     output_directory_path='../Ranked_Accuracy/HbCD')

    print("Creating submission collection for HbCD, including non-ranked")
    collection_hbcd_nonranked = HostGuestSubmissionCollection(submissions_hbcd, experimental_data,
            output_directory_path='../All_Accuracy/HbCD',
            ignore_refcalcs = False, ranked_only = False)
    
    # Create ranked submission for combine set of hosts. Will be for ranked molecule statistics
    print("Creating submission collection for combined set of hosts")
    collection_all = HostGuestSubmissionCollection(submissions_all, experimental_data,
                                                   output_directory_path='../Ranked_Accuracy/MoleculesStatistics', allow_multiple = True)

    # Create ranked and non-ranked submission collection for combine set of hosts (all submissions). Will be for molecule statistics of all methods. 
    #print("Creating submission collection for combined set of hosts, including non-ranked")
    collection_all_nonranked = HostGuestSubmissionCollection(submissions_all, experimental_data, 
            output_directory_path='../All_Accuracy/MoleculesStatistics', 
            allow_multiple = True, ignore_refcalcs = False, ranked_only = False)

    # Systems to be excluded (optionals and systems not detected).Currently only WP6-G4.
    def remove_optional(submission_collection_data):
        return submission_collection_data[(submission_collection_data.system_id != 'WP6-G4')]

    #make new collections and remove optionals. Currently only applies to WP6-G4.
    print("Making new collection set (ranked only) & removing optional host-guest systems from WP6 collection")
    collection_wp6_no_optional = HostGuestSubmissionCollection(submissions_wp6, experimental_data,
                                                                    output_directory_path='../Ranked_Accuracy/WP6_no_optional')
    collection_wp6_no_optional.data = remove_optional(collection_wp6_no_optional.data)

    print("Making new collection set (including nonranked) & removing optional host-guest systems from WP6 collection")
    collection_wp6_nonranked_no_optional = HostGuestSubmissionCollection(submissions_wp6, experimental_data, 
            output_directory_path='../All_Accuracy/WP6_no_optional', 
            ignore_refcalcs = False, ranked_only = False)
    collection_wp6_nonranked_no_optional.data = remove_optional(collection_wp6_nonranked_no_optional.data)

    print("Removing optional host-guest systems from entire ranked collection")
    collection_all.data = remove_optional(collection_all.data)

    print("Removing optional host-guest systems from entire ranked+nonranked collection")
    collection_all_nonranked.data = remove_optional(collection_all_nonranked.data)


    # =============================================================================
    # CREATE AUTOMATIC ANALYSIS ON THE REPO.
    # =============================================================================

    sns.set_style('whitegrid')

    # NOTE: Do not include collection_all or collection_all_nonranked here. 
    # Generate correlation plots and statistics
    for collection in [collection_cd, collection_cd_nonranked, collection_bcd, collection_bcd_nonranked, collection_hbcd, collection_hbcd_nonranked, collection_wp6_no_optional, collection_wp6_nonranked_no_optional, collection_wp6, collection_wp6_nonranked]:
        sns.set_context('notebook')
        collection.generate_correlation_plots()

        sns.set_context('talk')
        caption = ''
        collection.generate_statistics_tables(stats_funcs, subdirectory_path='StatisticsTables',
                                              groupby='name', extra_fields=['sid'],
                                              sort_stat='RMSE', ordering_functions=ordering_functions,
                                              latex_header_conversions=latex_header_conversions,
                                              caption=caption)

        sns.set_context('paper', font_scale=0.7)
        collection.plot_bootstrap_distributions(stats_funcs, subdirectory_path='StatisticsPlots',
                                                groupby='name', ordering_functions=ordering_functions,
                                                latex_header_conversions=latex_header_conversions,
                                                stats_limits=stats_limits)

    # Generate molecule statistics and plots for all ranked submissions (for now only WP6) 
    # Don't modify original collection_all as we'll use it later.
    collection = copy.deepcopy(collection_all)
    
    # Exclude "OPTIONAL" HOST-GUEST SYSTEMS
    collection.data = collection.data[~collection.data.system_id.isin({'WP6-G4'})]
    collection.generate_molecules_plot()
    collection.generate_statistics_tables(stats_funcs_molecules, 'StatisticsTables', groupby='system_id',
                                          sort_stat='MAE', ordering_functions=ordering_functions,
                                          latex_header_conversions=latex_header_conversions)
    collection.plot_bootstrap_distributions(stats_funcs_molecules, subdirectory_path='StatisticsPlots',
                                            groupby='system_id', ordering_functions=ordering_functions,
                                            latex_header_conversions=latex_header_conversions)

    # Generate molecule statistics and plots for all submissions (ranked and non ranked)
    collection_nr = copy.deepcopy(collection_all_nonranked)
    # Exclude "OPTIONAL" HOST-GUEST SYSTEMS
    collection_nr.data = collection_nr.data[~collection_nr.data.system_id.isin({'WP6-G4'})]
    collection_nr.generate_molecules_plot()
    collection_nr.generate_statistics_tables(stats_funcs_molecules, 'StatisticsTables', groupby='system_id', 
            sort_stat='MAE', ordering_functions=ordering_functions, 
            latex_header_conversions=latex_header_conversions)
    collection_nr.plot_bootstrap_distributions(stats_funcs_molecules, subdirectory_path='StatisticsPlots', 
            groupby='system_id', ordering_functions=ordering_functions,
            latex_header_conversions=latex_header_conversions)



    # =============================================================================
    # FIGURES AND TABLES GENERATED FOR THE PAPER
    # =============================================================================

    # Regenerate submission collection, this time without discarding
    # the methods that were applied to only one of the two sets.
    #submissions_oa_temoa = load_submissions(HostGuestSubmission, HOST_GUEST_OA_SUBMISSIONS_DIR_PATH, user_map)
    #submissions_oa_temoa = merge_submissions(submissions_oa_temoa, discard_not_matched=False)
    #collection_oa_temoa = HostGuestSubmissionCollection(submissions_oa_temoa, experimental_data,
    #                                                    output_directory_path='../Accuracy/OA-TEMOA')

    # Create a set of all the methods.
    all_methods = set(collection_all_nonranked.data.method.unique())
    #all_methods.update(set(collection_cb8.data.method.unique()))
    #all_methods.update(set(collection_temoa.data.method.unique()))
    #all_methods.update(set(collection_temoa_nonranked.data.method.unique()))
    #all_methods.update(set(collection_teetoa.data.method.unique()))
    #all_methods.update(set(collection_teetoa_nonranked.data.method.unique()))

    # Create a set of only ranked methods
    ranked_methods = set(collection_all.data.method.unique())

    # Submissions using experimental corrections.
    #is_corrected = lambda m: ('MovTyp' in m and m[-1] != 'N') or 'MMGBSA' in m or 'MMPBSA' in m 
    #corrected_methods = {m for m in all_methods if is_corrected(m)}

    # For movable type we plot only GE3N, GE3O, GE3L, KT1N, KT1L, GT1N, GT1L.
    exclusions = {} 
    #exclusions = {
    #    'MovTyp-GD1N', 'MovTyp-GD1O', 'MovTyp-GD1L', 'MovTyp-GD3N', 'MovTyp-GD3L',
    #    'MovTyp-GE3S', 'MovTyp-GE3U', 'MovTyp-GE3Z', 'MovTyp-GT1O', 'MovTyp-GT3N',
    #    'MovTyp-GT3S', 'MovTyp-GT3U', 'MovTyp-GT3L', 'MovTyp-GU1N', 'MovTyp-GU1O',
    #    'MovTyp-GU1L', 'MovTyp-GU3N', 'MovTyp-GU3L'
    #}


    # FIGURE 2: Figure experiment distributions.
    # -------------------------------------------
    sns.set_style('darkgrid')
    sns.set_context('paper', font_scale=1.0)
    host_names = [system_id.split('-')[0] for system_id in experimental_data.index]
    system_id_idx = list(range(len(host_names)))
    data = experimental_data.assign(host_name=pd.Series(host_names).values,
                                    system=system_id_idx)
    data = data.rename(columns={
        '$\Delta$G': '$\Delta$G [kcal/mol]',
        'host_name': 'Host name',
        'system': 'Host-guest system'
    })
    fix, ax = plt.subplots(figsize=(4.5, 4.5))
    ax = sns.stripplot(y=data.index, x='$\Delta$G [kcal/mol]', data=data,
                       hue='Host name', size=9, palette=HOST_PALETTE, ax=ax)
    ax.errorbar(y=list(range(len(data.index))), x=data['$\Delta$G [kcal/mol]'].values,
                xerr=data['d$\Delta$G'].values, fmt='none', elinewidth=1, ecolor='black',
                capsize=2, capthick=1, zorder=10)
    ax.set_xlim((-15, 0))
    plt.tight_layout(pad=0.2)
    # plt.show()
    plt.savefig('../Ranked_Accuracy/PaperImages/Figure2_experimental_measurements.pdf')


    # FIGURE 3: Figure correlation plots free energies.
    # --------------------------------------------------
    sns.set_style('whitegrid')
    sns.set_context('paper', font_scale=0.8)

    def correlation_plots(plotted_methods, file_name):
        """Shortcut to create correlation plots."""
        n_methods = len(plotted_methods)

        if n_methods > 4:
            n_cols = 5
            n_rows = int(np.floor(n_methods/(n_cols-1)))
        else:
            n_cols = n_methods
            n_rows = 1
        plot_size = 15.25 / n_cols
        fig = plt.figure(figsize=(n_cols*plot_size, n_rows*plot_size))
        grid = plt.GridSpec(nrows=n_rows, ncols=n_cols*3)
        # All rows have 5 plots.
        axes = []
        #TODO: This needs generalization
        for row_idx in range(n_rows-1):
            axes.extend([fig.add_subplot(grid[row_idx, c:c+2]) for c in range(1,10,2)])
            #axes.extend([fig.add_subplot(grid[row_idx, c:c+2]) for c in range(1,7,2)])
        axes.extend([fig.add_subplot(grid[-1, c:c+2]) for c in range(1,5,2)])

        # Associate a color to each host.
        for method, ax in zip(plotted_methods, axes):
            # Isolate statistics of the method.
            data = collection_all.data[collection_all.data.method == method]
            # Build palette.
            palette = [HOST_PALETTE[host_name] for host_name in sorted(data.host_name.unique())]
            # Add color for regression line over all data points.
            palette += [HOST_PALETTE['other1']]
            # Plot correlations.
            plot_correlation(x='$\Delta$G (expt) [kcal/mol]', y='$\Delta$G (calc) [kcal/mol]',
                             data=data, title=method, hue='host_name', color=palette,
                             shaded_area_color=HOST_PALETTE['other2'], ax=ax)
            # Remove legend and axes labels.
            ax.legend_.remove()
            ax.set_xlabel('')
            ax.set_ylabel('')
            # Make title and axes labels closer to axes.
            #ax.set_title(ax.get_title(), pad=1.5)
            ax.set_title(ax.get_title(), pad=2.0)
            ax.tick_params(pad=3.0)

        # Use a single label for the figure.
        fig.text(0.015, 0.5, '$\Delta$G (calc) [kcal/mol]', va='center', rotation='vertical', size='large')
        fig.text(0.5, 0.015, '$\Delta$G (exp) [kcal/mol]', ha='center', size='large')

        plt.tight_layout(pad=0.9, rect=[0.0, 0.025, 1.0, 1.0])
        #plt.tight_layout(pad=1.5)
        plt.savefig('../Ranked_Accuracy/PaperImages/{}.pdf'.format(file_name))
    
    # Create correlation plot of ranked methods only. 
    correlation_plots(
        #plotted_methods = sorted(set(all_methods) - set(exclusions) - {'NULL'}),
        #plotted_methods = sorted(set(all_methods) - set(exclusions)),
        plotted_methods = sorted(set(ranked_methods)),
        file_name='Figure_correlation_plots_ranked_methods'
    )

    # Supplementary figure with correlations plots of both ranked and non-ranked methods
    #correlation_plots(
    #    plotted_methods = sorted(set(all_methods)),
    #    file_name='Figure_correlation_plots_all_methods'
    #)


    # FIGURE 5: Figure statistics by molecule.
    # -----------------------------------------
    sns.set_style('whitegrid')
    sns.set_context('paper', font_scale=1.0)

    def get_errs(stat_name, data):
        stat_name_lb = stat_name + '_lower_bound'
        stat_name_ub = stat_name + '_upper_bound'
        errs = []
        for i, x in data[[stat_name, stat_name + '_lower_bound', stat_name + '_upper_bound']].iterrows():
            errs.append([abs(x[stat_name_lb] - x[stat_name]), abs(x[stat_name_ub] - x[stat_name])])
        return np.array(list(zip(*errs)))

    stats_names = ['RMSE', 'ME']
    fig, axes = plt.subplots(ncols=len(stats_names), figsize=(7.25/3.1*2, 5), sharey=True)
    # only ranked methods
    statistics = pd.read_json('../Ranked_Accuracy/MoleculesStatistics/StatisticsTables/statistics.json', orient='index')

    # Remove OPTIONALS (). Plots error by molecule only for ranked
    #statistics = statistics[~statistics.index.isin({'WP6-G4'})]
    statistics.sort_values(by='RMSE', inplace=True)
    for ax, stats_name in zip(axes, stats_names):
        # Build palette.
        palette = [HOST_PALETTE[guest_name.split('-')[0]] for guest_name in statistics.index.values]
        # Convert errors for printing.
        rmse_errs = get_errs(stats_name, statistics)
        ax = sns.barplot(x=stats_name, y=statistics.index, data=statistics, xerr=rmse_errs,
                         palette=palette, ax=ax)
        if stats_name == 'ME':
            ax.set_xlim((-3, 8))
        ax.set_xlabel('$\Delta$G ' + stats_name + ' [kcal/mol]')

    plt.tight_layout(pad=0.3)
    # plt.show()
    plt.savefig('../Ranked_Accuracy/PaperImages/error_by_molecule.pdf')

    # all methods including nonranked
    statistics_nr = pd.read_json('../All_Accuracy/MoleculesStatistics/StatisticsTables/statistics.json', orient='index')

    # Remove OPTIONALS (). Plots error by molecule for all submissions (includes non ranked)
    #statistics_nr = statistics_nr[~statistics_nr.index.isin({'WP6-G4'})]
    statistics_nr.sort_values(by='RMSE', inplace=True)
    for ax, stats_name in zip(axes, stats_names):
        # Build palette.
        palette = [HOST_PALETTE[guest_name.split('-')[0]] for guest_name in statistics_nr.index.values]
        # Convert errors for printing
        rmse_errs = get_errs(stats_name, statistics_nr)
        ax = sns.barplot(x=stats_name, y=statistics_nr.index, data=statistics_nr, xerr=rmse_errs,
                palette=palette, ax=ax)
        if stats_name == 'ME':
            ax.set_xlim((-3, 8))
        ax.set_xlabel('$\Delta$G ' + stats_name + ' [kcal/mol]')

    plt.tight_layout(pad=0.3)
    #plt.show()
    plt.savefig('../All_Accuracy/PaperImages/error_by_molecule.pdf')

    #Break before making next image. Test all of the above.
    import sys
    sys.exit('Stop before generating tightest binders plot')


    # Create table presenting which methods got the tightest binders.
    # ---------------------------------------------------------------

    sns.set_style('white')

    # Consider only the methods that are part of the main analysis.
    plotted_methods = sorted(set(ranked_methods) - set(exclusions))

    # Create a Dataframe summarizing if a method got the tightest binders correctly for each guest set.
    data = collections.OrderedDict()  # Data in dict format.
    # Tightest binders and columns of the Pandas Dataframe.
    tighest_binders =['WP6-G12']
    for method in plotted_methods:
        data[method] = []

        for collection, tighest_binder in zip([collection_cb8, collection_temoa, collection_teetoa], tighest_binders):
            method_data = collection.data[collection.data.method == method]

            # Get the tightest binder predicted by the method for the set.
            is_prediction_correct = 0
            is_prediction_missing = 0
            try:
                predicted_tightest_binder_id = method_data['$\Delta$G (calc) [kcal/mol]'].idxmin()
            except ValueError:
                # The method did not submit predictions for this guest set.
                is_prediction_missing = 1
            else:
                predicted_tightest_binder = method_data[method_data.index == predicted_tightest_binder_id]
                predicted_tightest_binder = predicted_tightest_binder.system_id.values[0]
                is_prediction_correct = 1 if predicted_tightest_binder == tighest_binder else 0

            # Push column in the table.
            data[method].append(is_prediction_correct)
            data[method].append(1 - is_prediction_correct - is_prediction_missing)
            data[method].append(is_prediction_missing)

    # Reorder records to display the methods that got the highest number of tighest binders right.
    def rank_method(d):
        # Correct predictions fully count.
        score = -sum(d[1][i] for i in range(0, len(d[1])) if i % 3 == 0)
        # Not submitted predictions count less.
        score -= sum(d[1][i]*0.1 for i in range(0, len(d[1])) if i % 3 == 2)
        # Secondary key is alphabetical.
        return (score, d[0])
    # Rank methods by
    data = sorted(data.items(), key=lambda d: rank_method(d), reverse=True)
    data = collections.OrderedDict(data)

    # Convert table to dataframe. 
    columns = ['WP6-G12', 'WP6-G12-incorrect', 'WP6-G12-notsubmitted']
    palette = [sns.desaturate(color, 0.75) for color in [HOST_PALETTE['WP6'], '0.7', 'white']]
    #palette = [sns.desaturate(color, 0.75) for color in [HOST_PALETTE['CB8'], '0.7', 'white',
    #                                                     HOST_PALETTE['TEMOA'], '0.7', 'white',
    #                                                     HOST_PALETTE['TEETOA'], '0.7', 'white']] 
    
    data = pd.DataFrame.from_dict(data, orient='index', columns=columns)

    # Plot percentage of correct binders across methods.
    for tighest_binder in tighest_binders:
        tot_predictions = len(data[tighest_binder]) - sum(data[tighest_binder + '-notsubmitted'])
        correct_predictions = sum(data[tighest_binder])
        print('{}: {}/{} ({:.2f}%)'.format(tighest_binder, correct_predictions, tot_predictions,
                                        correct_predictions/tot_predictions*100))

    # Plot table.
    #ax = data.plot.barh(stacked=True, color=palette, figsize=(7.25/3.1,5))
    ax = data.plot.barh(stacked=True, color=palette, figsize=(10, 10)) # TRY NEW FIG SIZE
    ax.xaxis.set_ticks([0.5, 1.5, 2.5, 3.5]) # ADD ANOTHER TICK 3.5
    ax.xaxis.set_ticklabels(['WP6']) # ADD ADDITIONAL HOSTS HERE
    ax.set_title('Methods predicting the tightest binders')

    # Configure lengend to hide the missing/incorrect labels.
    handles, labels = ax.get_legend_handles_labels()
    for i in reversed(range(len(handles))):
        if (('incorrect' in labels[i] and i > 2) or  # Leave only 1 'incorrect' label.
                ('notsubmitted' in labels[i])):  # Delete all 'incorrect' legend labels.
            del handles[i]
            del labels[i]
    # Change the unique label for 'notsubmitted' and move it at the end of the legend.
    labels[1] = 'Incorrect'
    handles.append(handles.pop(1))
    labels.append(labels.pop(1))
    # Locate legend on top of the plot.
    ax.legend(handles, labels, ncol=2, loc='upper right', markerfirst=False,
              bbox_to_anchor=(1.0, -0.03), columnspacing=1.5)

    plt.tight_layout(rect=[0, 0.0, 1, 1], pad=0.1)
    # plt.show()
    plt.savefig('../Ranked_Accuracy/PaperImages/tightest_binders.pdf')


    #Break before making next image. Test all of the above.
    #import sys
    #sys.exit('Stop before generating next plot/table')


    # Generate the initial statistics table for the paper.
    # -------------------------------------------------------
    stats_funcs = collections.OrderedDict([
         ('RMSE', rmse),
         ('ME', me),
         ('R2', r2),
         ('kendall_tau', kendall_tau)
     ])
    collection.generate_paper_table(stats_funcs, exclusions)

    #Break before making next image. Test all of the above
    #import sys
    #sys.exit('Stop before generating next plot/table')

    # Create initial table of methods with and without bonus challenge.
    # ------------------------------------------------------------------
    # We include only the ones that submitted entries for the bonus challenge.
    #exclusions = all_methods - {'ForceMatch', 'ForceMatch-QMMM', 'Tinker-AMOEBA',
    #                            'DFT(B3PW91)', 'DFT(B3PW91)-D3', 'MMPBSA-GAFF'}
    #collection = SplitBootstrapSubmissionCollection(collection_cb, collection_cb_no_bonus,
    #                                                hue='dataset', collection1_hue='BONUS', collection2_hue='NOBONUS',
    #                                                output_directory_path='../MergedCB8')
    #collection.generate_paper_table(stats_funcs, exclusions)

    # # Compute statistics similar to old review papers.
    # # ----------------------------------------
    # def print_stats(methods, stats_names, absolute_statistics, relative_statistics):
    #     for method in methods:
    #         print('{:15s} & {:13s} '.format(method, 'OA/TEMOA'), end='')
    #         for statistics, suffix in [(absolute_statistics, ''), (relative_statistics, '_o')]:
    #             for i, stats_name in enumerate(stats_names):
    #                 stats_name += suffix
    #                 stats, ci, bootstrap_distribution = statistics[method][stats_name]
    #                 # print('{}-{}: {:.1f} ({:.1f} +- {:.1f})'.format(method, stats_name, stats,
    #                 #                                                 np.mean(bootstrap_distribution),
    #                 #                                                 np.std(bootstrap_distribution)))
    #                 # end = ' \\\\\n' if (i == len(stats_names) - 1 and suffix == '_o') else ' & '
    #                 print(' & {:.1f} ({:.1f} $\pm$ {:.1f})'.format(stats, np.mean(bootstrap_distribution),
    #                                                          np.std(bootstrap_distribution)), end='')
    #         print('     \\\\')
    #
    # # Compute offset statistics
    # import scipy
    #
    # def rmse_offset(data):
    #     mse = me(data)
    #     x, y = data.T
    #     offset_error = np.array(x) - np.array(y) - mse
    #     rmse = np.sqrt((offset_error**2).mean())
    #     return rmse
    #
    # def kendall_tau_offset(data):
    #     mse = me(data)
    #     x, y = data.T
    #     y += mse
    #     correlation, p_value = scipy.stats.kendalltau(x, y)
    #     return correlation
    #
    # def r2_offset(data):
    #     mse = me(data)
    #     x, y = data.T
    #     y += mse
    #     slope, intercept, r_value, p_value, stderr = scipy.stats.linregress(x, y)
    #     return r_value**2
    #
    # stats_funcs_offset = collections.OrderedDict([
    #     ('RMSE_o', rmse_offset),
    #     ('R2_o', r2_offset),
    #     ('kendall_tau_o', kendall_tau_offset)
    # ])
    #
    # offset_methods = ['DFT(TPSS)-D3', 'MovTyp-KT1N', 'MovTyp-KT1L', 'MovTyp-GE3N',
    #                   'SOMD-A-nobuffer', 'SOMD-C-nobuffer','SOMD-D-nobuffer']
    # collection = copy.deepcopy(collection_oa_temoa)
    # # offset_methods = ['Tinker-AMOEBA']
    # # collection = copy.deepcopy(collection_cb_no_bonus)
    #
    # # Obtain absolute statistics
    # stats_names, stats_funcs = zip(*stats_funcs.items())
    # cache_file_path = os.path.join(collection.output_directory_path, 'bootstrap_distributions.p')
    # absolute_statistics = collection._get_bootstrap_statistics('method', stats_names,
    #                                                            stats_funcs, cache_file_path)
    #
    # # Discard methods that we don't need to compute the offset statistics of and obtain offset statistics.
    # collection.data = collection.data[collection.data.method.isin(offset_methods)]
    # stats_names_offset, stats_funcs_offset = zip(*stats_funcs_offset.items())
    # cache_file_path = os.path.join(collection.output_directory_path, 'bootstrap_distributions_offset.p')
    # offset_statistics = collection._get_bootstrap_statistics('method', stats_names_offset,
    #                                                          stats_funcs_offset, cache_file_path)
    # print_stats(methods=offset_methods, stats_names=['RMSE', 'R2', 'kendall_tau'],
    #             absolute_statistics=absolute_statistics, relative_statistics=offset_statistics)
