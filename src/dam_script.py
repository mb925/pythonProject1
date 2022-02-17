# Dependencies
import pandas as pd
import numpy as np
import time
import json
import sys
import re


# Define summary table
class Summary(object):
    """ Summary table for MobiDB

    The summary table for MobiDB is stored as a list of dictionaries within a JSON file.
    Each entry is associated to a unique UniProt accession. It contains attributes of different types: some are descriptive,
    some others map information between external resources and the UniProt sequence. Such mappings represent features.

    Features are in turn identified by the following triplet: evidence, feature and source. Evidence defines the quality of
    the annotation (e.g. curated, derived, homology, prediction). Feature defines the type of annotation (e.g. disorder, lip, ...).
    Source defines the annotation source which can either be a software or an external database (e.g. mobidb_lite, disprot, merged,
    th_90, ...).

    The evidence field can assume the following values: `curated` is the most reliable one and defines curated experimental annotations
    from partner databases; `derived` is less reliable and determines entries automatically derived from primary data; `homology` is
    more comprehensive and identifies annotations propagated by aligning curated regions available from Ensembl and Gene Trees;
    `prediction` is the most comprehensive and defines annotations retrieved by running software tools agains sequences in UniProt.

    Various consensus levels are available for each mapping between UniProt sequences and their features. `th_90` defines regions for
    which there is an agreement of 90% at least; `th_50` is less strict and defines regions where the mjority of annotations agree;
    `merge` is the most permessive since it just merges the different annotations; `priority` combines together annotations from
    different sources and with different levels in a hierarchical manner (e.g. `curated` overwrites `derived`, ...).

    The main features that can be found in each entry of the MobiDB summary are about intrinsically disordered regions (IDR) and
    linearly interacting peptides (LIP). LIP are regions within IDR that can perform interactions through short peptide motifs.

    There are several sources for IDR, either manually curated or automatically annotated. The most important ones are listed below:
    - `curated-disorder-<db>` defines disorder identified on a manually curated database (e.g. `disprot`, `ideal`, 'uniprot`, ...);
    - `homology-disorder-<db>` defines disorder identified by homology searching on a database (e.g. `disprot`, `ideal`, `uniprot`);
    - `prediction-disorder-mobidb_lite` identifies regions preditced to be IDR by the MobiDB-lite software;
    - `prediction-<subregion>-mobidb_lite_sub` identifies specific sub-regions within a region prodicted to be an IDR by the
    MobiDB-lite software (e.g. `prediction-low_complexity-mobidb_lite_sub`, `prediction-cystein_rich-mobidb_lite_sub`, ...).

    Among the sources available for linear interacting peptides (LIP) there are the following:
    - `curated-lip-<db>` identifies manually curated LIP found in other sources (e.g. `curated-lip-disprot`, `curated-lip-dibs`);
    - `derived-lip-<consensus>` identifies regions predicted to be LIP by the Fast Linear Interacting Peptides Predictor (FLIPPER)
    software. Such predictions are derived from structure and hence are thought to be sufficiently reliable;
    - `prediction-lip-anchor` identifies regions predicted to be LIP by the anchor software. It uses only sequence information and
    is not considered as reliable as FLIPPER.

    Attributes
    ----------
    entries : dict
        A dictionary of Summary.Entry classes whose keys are Entry.accession attributes;

    Methods
    -------
    load_json : Summary
        Load JSON entries and features from JSON file
    parse_json : Summary
        Read entries and features from JSON formatted list
    get_entries : pd.DataFrame
        Returns the DataFrame of entries
    get_features : pd.DataFrame
        Returns the DataFrame of features
    """

    # Constructor
    def __init__(self, entries):
        # Set attributes
        self.entries = entries

    def __len__(self):
        return len(self.entries)

    @classmethod
    def load_json(cls, file_or_path):
        # Check if input handler is path
        is_path = isinstance(file_or_path, str)
        # Initialize file handler
        handler = open(file_or_path, 'r') if is_path else file_or_path
        # Load entries from JSON
        count, summary = cls.parse_json(json.load(handler))
        # If opened here: close connection to file
        if not is_path: handler.close()
        # Return parsed summary
        return count, summary

    @classmethod
    def load_ndjson(cls, file_or_path):
        # Check if input handler is path
        is_path = isinstance(file_or_path, str)
        # Initialize file handler
        handler = open(file_or_path, 'r') if is_path else file_or_path
        # Load entries from newline-delimited JSON (ndjson)
        count, summary = cls.parse_json(map(json.loads, handler))
        # If opened here: close connection to file
        if not is_path: handler.close()
        # Return parsed summary
        return count, summary

    @classmethod
    def parse_json(cls, entries):
        # Initialize entries counter
        num_entries = 0
        # Initialize entries dictionary
        entries_dict = dict()
        # Loop through each entry in input dictionary
        for _, entry in enumerate(entries):
            # Parse entry from JSON
            entry = cls.Entry.parse_json(entry)
            # Store entry
            entries_dict.setdefault(entry.acc, entry)
            # Update counter
            num_entries += 1
        # Return counter and summary object
        return num_entries, Summary(entries_dict)

    # Cast entries to dataframe
    def to_dataframe(self):
        # Initialize list of identities for UniProt sequences
        uniprot_sequences = list()
        # Loop through each entry
        for _, curr_entry in self.entries.items():
            # Store current entry
            uniprot_sequences.append(curr_entry.to_dict())
        # Make main table for UniProt sequences
        return pd.DataFrame(uniprot_sequences)

    class Entry(object):

        # Define allowed features by regular expression
        re_features = [r'derived-lip-*',  # Allow all curated LIP predictions from FLIPPER
                       r'curated-disorder-*',  # Allow all curated disorder predictions
                       r'homology-domain-*']  # Allow domains detected by homology

        def __init__(self, accession, sequence, length, uniref100='', uniref90='', uniref50='', uniparc='', name='',
                     organism='', reviewed=False, features=dict()):
            # Store attributes
            self.accession = accession
            self.sequence = sequence
            self.length = length
            self.uniref100 = uniref100
            self.uniref90 = uniref90
            self.uniref50 = uniref50
            self.uniparc = uniparc
            self.name = name
            self.organism = organism
            self.reviewed = reviewed
            self.features = features

        @property
        def acc(self):
            return self.accession

        @property
        def seq(self):
            return self.sequence

        @property
        def feature(self):
            # Define a feature describing the whole sequence
            return Summary.Feature(
                # Reference itself
                entry=self,
                # Copy some attributes
                accession=self.accession,
                name=self.name,
                # Set sequence length as the whole sequence
                begin=1, end=len(self) + 1
            )

        # Counts letters in sequence string
        def __len__(self):
            return len(self.sequence)

        @classmethod
        def parse_json(cls, entry):
            # Parse attributes
            accession = str(entry.get('acc', None))
            sequence = str(entry.get('sequence', ''))
            length = int(entry.get('length', 0))
            uniref100, uniref90, uniref50 = tuple(str(entry.get('uniref%d' % i, '')) for i in (100, 90, 50))
            uniparc = str(entry.get('uniparc', ''))
            name = str(entry.get('name', ''))
            organism = str(entry.get('organism', ''))
            reviewed = bool(entry.get('reviewed', False))
            # Initialize features dictionary
            features = dict()
            # Generate parsed entry
            parsed = cls(accession, sequence, length, uniref100, uniref90, uniref50, uniparc, name, organism, reviewed,
                         features)
            # Loop through each feature
            for key, prediction in entry.items():
                # Initialize allowed regex flag
                is_allowed = False
                # Split key into pieces using dash separator
                pieces = re.split(r'-', key)
                # Loop through each allowed regular expression regex
                for regex in cls.re_features:
                    # Check regex
                    is_allowed = is_allowed or bool(re.match(regex, key))
                    # # Debug
                    # print('Debug', '%s - %s - %s' % (key, regex, str(is_allowed)))
                # Skip if feature is not allowed
                if not is_allowed:
                    continue
                # Split key into (evidence, name, source, ...) tuple
                feature_id = tuple(pieces)
                # Initialize list for current feature
                features_list = features.setdefault(feature_id, list())
                # Get regions
                regions = prediction.get('regions', list())
                regions_ids = prediction.get('regions_ids', list())
                regions_names = prediction.get('regions_names', list())
                # # If region names of IDs are not specified, skip them
                # if 'regions_ids' not in prediction.keys():
                #     continue
                # elif 'regions_names' not in prediction.keys():
                #     continue
                # Parse each region
                for i, _ in enumerate(regions):
                    # Get region bounderies on entry's
                    region_begin, region_end = tuple(regions[i])
                    # Get region id
                    region_id = regions_ids[i] if (i < len(regions_ids)) else ''
                    # Get region name
                    region_name = regions_names[i] if (i < len(regions_names)) else ''
                    # Create new feature for the region itself
                    feature = Summary.Feature(accession=region_id, name=region_name,
                                              begin=region_begin, end=region_end,
                                              entry=parsed)
                    # Update features list
                    features_list.append(feature)
            # Return new entry
            return parsed

        def to_dict(self):
            # Cast entry to dictionary
            return {'entry_acc': self.acc,
                    'entry_seq': self.seq,
                    'entry_len': len(self),
                    'entry_name': self.name,
                    'entry_organism': self.organism,
                    'entry_uniref100': self.uniref100,
                    'entry_uniref90': self.uniref90,
                    'entry_uniref50': self.uniref50,
                    'entry_uniparc': self.uniparc,
                    'entry_reviewed': self.reviewed}

        @property
        def num_lips(self):
            # Just call feature's method on wholse sequence
            return self.feature.num_lips

        @property
        def num_disorder(self):
            # Just call feature's method on whole sequence
            return self.feature.num_disorder

        @property
        def num_pfam(self):
            # Just call feature's method on whole sequence
            return self.feature.num_pfam

        # def get_derived_lip(self, *args, **kwargs):
        #     return self.feature.get_derived_lips(*args, **kwargs)

        def get_link_domains(self, *args, **kwargs):
            return self.feature.get_link_domains(*args, **kwargs)

    class Feature(object):

        # (Static) evidence levels
        EVIDENCE_CURATED = 'curated'
        EVIDENCE_DERIVED = 'derived'
        EVIDENCE_HOMOLOGY = 'homology'
        EVIDENCE_PREDICTION = 'prediction'

        # (Static) consensus
        CONSENSUS_STRICT = 'th_90'
        CONSENSUS_MAJORITY = 'th_50'
        CONSENSUS_MERGE = 'merge'
        CONSENSUS_PRIOROTY = 'priority'

        # (Static) names
        NAME_LIP = 'lip'
        NAME_PFAM = 'pfam'
        NAME_DISORDER = 'disorder'

        # (Static) mapping between source name and code
        feature_sources = {'strict': CONSENSUS_STRICT,
                           'majority': CONSENSUS_MAJORITY,
                           'merge': CONSENSUS_MERGE,
                           'priority': CONSENSUS_PRIOROTY}

        def __init__(self, entry, accession='', name='', begin=1, end=0):
            # Store attributes
            self.entry = entry  # Reference to entry where domain was found
            self.accession = accession  # Identifier of the region
            self.name = name  # Name of the region
            self.begin = begin  # Index on the entry's sequence where the domain begins
            self.end = end  # Index on the entry's sequence where the domain ends

        @property
        def features(self):
            return self.entry.features

        @property
        def i(self):
            # Returns actual index starting from 0
            return self.begin - 1

        @property
        def j(self):
            return self.end

        @property
        def acc(self):
            return self.accession

        def __len__(self):
            return self.j - self.i

        @property
        def bin_sequence(self):
            # Get shape of current entry
            m = len(self.entry)
            # Initialize binary sequence
            bin_sequence = np.zeros((m,), dtype=int)
            # Update binary sequence
            bin_sequence[self.i:self.j] = 1
            # Return binary sequence
            return bin_sequence

        @property
        def num_lips(self):
            # Initialize number of LIPs for various sources
            num_lips = dict()
            # Retrieve sources for features
            feature_sources = self.feature_sources
            # Loop through allowed features
            for source_name, source_code in feature_sources.items():
                # Define feature's identifier
                feature_id = self.EVIDENCE_DERIVED, self.NAME_LIP, source_code
                # Get features
                features_list = self.features.get(feature_id, list())
                # Get LIP consensus
                lip_consensus = self.get_consensus(*features_list)
                # Compute number of LIP residues for current source
                num_lips.setdefault(source_name, lip_consensus.sum())
            # Return number of LIPs for each consensus level
            return num_lips

        @property
        def num_disorder(self):
            """ Get number of disordered residues over current region

            Returns
            -------
            num_disorder : dict
                Dictionary mapping consensus level to number of disordered residues on current region
            """
            # Initialize dictionary mapping consensus to disorder statistics
            num_disorder = dict()
            # Retrieve feature sources
            feature_sources = {**self.feature_sources, **{'ideal': 'ideal', 'disprot': 'disprot'}}
            # Loop through each available source
            for source_name, source_code in feature_sources.items():
                # Define current feature's identifier
                feature_id = self.EVIDENCE_CURATED, self.NAME_DISORDER, source_code
                # Get features
                features_list = self.features.get(feature_id, list())
                # Compute disorder consensus
                disorder_consensus = self.get_consensus(*features_list)
                # Update number of disordered proteins for current source
                num_disorder.setdefault(source_name, disorder_consensus.sum())
            # Return dictionary mapping consensus to disorder statistics
            return num_disorder

        @property
        def num_pfam(self):
            # Get shape of current entry
            m = len(self.entry)
            # Initialize domain sequence
            sequence_bin = np.zeros((m,), dtype=int)
            # Define feature identifier
            feature_id = 'homology', 'domain', 'pfam'
            # Get list of Pfam domains
            pfam_domains = self.features.get(feature_id, list())
            # Initialize number of pfam residues and domains
            num_pfam_domains, num_pfam_residues = 0, 0
            # Loop through each LIP region
            for pfam_region in pfam_domains:
                # Skip if region is left
                if pfam_region.j <= self.i:
                    continue
                # Skip if region is right
                if pfam_region.i >= self.j:
                    continue
                # Update number of pfam domains
                num_pfam_domains += 1
                # Define LIP regions on binary sequence
                sequence_bin[pfam_region.i:pfam_region.j] = 1
            # Get binary sub-sequence
            feature_bin = sequence_bin[self.i:self.j]
            # Compute number of Pfam residues for current entry
            num_pfam_residues = feature_bin.sum()
            # Return both number of pfam domains and pfam residues
            return num_pfam_domains, num_pfam_residues

        def to_dict(self):
            # Cast feature to dictionary
            return {'entry_acc': self.entry.acc,
                    'region_acc': self.acc,
                    'region_name': self.name,
                    'region_begin': self.begin,
                    'region_end': self.end,
                    'region_len': len(self)}

        # Retrieve coverage on another feature (same entry)
        def get_coverage(self, feature, on='self'):
            # check whether input feature is of correct type
            if not isinstance(feature, Summary.Feature):
                raise ValueError('Only Feature instances are accepted')
            # Check whether features are on the same entry
            if self.entry != feature.entry:
                raise ValueError('Only features on the same sequence can be compared')
            # Get binary sequences
            bin_feature = feature.bin_sequence
            bin_sequence = self.bin_sequence
            # Compute overlap between sequences
            bin_overlap = bin_sequence * bin_feature
            # Initialize numerator and denominator
            num, den = .0, .0
            # Case coverage of the other feature on itself
            if on == 'self':
                num = np.sum(bin_overlap[bin_sequence > 0])
                den = np.sum(bin_sequence)
            # Case coverage of itself on another feature
            elif on == 'other':
                num = np.sum(bin_overlap[bin_feature > 0])
                den = np.sum(bin_feature)
            # Case coverage of both features together
            elif on == 'both':
                num = np.sum(bin_overlap)
                den = np.sum((bin_sequence + bin_feature) > 0)
            # Otherwise, no valid method has been selected
            else:
                raise ValueError('No valid coverage method has been selected')
            # Return coverage
            return num / den

        # Retrieve consensus of one or more other features
        def get_consensus(self, *features):
            # Define length of sequence
            m = len(self.entry)
            # Initialize empty binary sequence
            bin_sequence = np.zeros((m,), dtype=int)
            # Loop through each given feature in array
            for curr_feature in features:
                # Check if item is a feature
                if not isinstance(curr_feature, Summary.Feature):
                    raise ValueError('Only Feature instances are accepted')
                # Check if entry is the same
                if self.entry is not curr_feature.entry:
                    raise ValueError('Only features on the same sequence can be compared')
                # Retrieve annotation boundaries
                i, j = curr_feature.i, curr_feature.j
                # Update binary sequence
                bin_sequence[i:j] += 1
            # Retrieve consensus on current feature
            return bin_sequence[self.i:self.j]

        def get_link_domains(self, *args, **kwargs):
            """ Retrieve regions bound by Pfam domains

            Returns
            -------
            link_domains : list
                List of regions on UniProt sequence which are bound by Pfam domains.
                The accession field of each of such region is the union between left
                and right bounding domains.
            """
            # Initialize list of linking domains
            link_domains = list()
            # Initialize key of feature holding Pfam domains on current UniProt sequence
            feature_id = 'homology', 'domain', 'pfam'
            # Get list of Pfam domain regions on current UniProt sequence
            pfam_domains = self.features.get(feature_id, list())
            # Sort domains by starting index
            pfam_domains.sort(key=lambda r: r.i, reverse=False)
            # Loop from second to last Pfam domain
            for k in range(1, len(pfam_domains)):
                # Get left and right bounding domains
                left_domain, right_domain = pfam_domains[k - 1], pfam_domains[k]
                # Get left and right bounding indices
                i, j = left_domain.j, right_domain.i
                # Skip if bounded domain is too small
                if j <= i:
                    continue
                # Skip if bounded domain is out of current domain
                if (j <= self.i) or (self.j <= i):
                    continue
                # Otherwise, define a boudaries for linking domain
                begin, end = (left_domain.end + 1), (right_domain.begin - 1)
                # Join together UniProt identifier and Pfam identifiers
                accession = '%s_%s' % (left_domain.acc, right_domain.acc)
                # Then, instantiate linking domain region
                curr_domain = self.__class__(
                    entry=self.entry, begin=begin, end=end, accession=accession,
                    name='Linking domain between domains %s and %s' % (left_domain.acc, right_domain.acc)
                )
                # Store linking domain
                link_domains.append(curr_domain)
            # Return list of linking domains
            return link_domains

        # TODO Get features overlapping the current one
        def get_inner_regions(self, feature_sources=dict()):
            raise NotImplementedError

        def get_lip_regions(self, *args, **kwargs):
            """ Retrieve LIP regions within current feature

            Returns
            -------
            regions_dict : dict
                Dictionary mapping consensus levels to list of features. Each list
                contains as many features as the number of LIP regions with given
                consensus level.
            """
            # Initialize regins dictionary
            regions_dict = dict()
            # Define current entry
            entry = self.entry
            # Define feature sources
            feature_sources = self.feature_sources
            # Loop through each feature source
            for feature_name, consensus_level in feature_sources.items():
                # Define feature identifier
                feature_id = 'derived', 'lip', consensus_level
                # Initialize list for current consensus level
                regions_list = regions_dict.setdefault(feature_name, list())
                # Loop through each LIP region for current consensus level
                for curr_region in self.features.get(feature_id, list()):
                    # Define accession for current region
                    accession = curr_region.acc
                    # Define name for current region
                    name = curr_region.name
                    # Define boundaries of current region
                    begin = max(curr_region.begin, self.begin)
                    end = min(curr_region.end, self.end)
                    # Skip if current region does not overlap
                    if end < begin:
                        continue
                    # Otherwise, increase the list for current consensus level
                    regions_list.append(self.__class__(
                        accession=accession, name=name,
                        entry=entry, begin=begin, end=end
                    ))
            # Define regions dictionary
            return regions_dict

        def get_disordered_regions(self, *args, **kwargs):
            """ Retrieve disordered regions within the current feature

            Returns
            -------
            regions_dict : dict
                Dictionary mapping consensus levels to list of features. Each list
                contains as many features as the number of disordered regions with given
                consensus level.
            """
            # Initialize regins dictionary
            regions_dict = dict()
            # Define current entry
            entry = self.entry
            # Define feature sources
            feature_sources = {**self.feature_sources, **{'ideal': 'ideal', 'disprot': 'disprot'}}
            # Loop through each feature source
            for feature_name, consensus_level in feature_sources.items():
                # Define feature identifier
                feature_id = self.EVIDENCE_CURATED, self.NAME_DISORDER, consensus_level
                # Initialize list for current consensus level
                regions_list = regions_dict.setdefault(feature_name, list())
                # Loop through each LIP region for current consensus level
                for curr_region in self.features.get(feature_id, list()):
                    # Define accession for current region
                    accession = curr_region.acc
                    # Define name for current region
                    name = curr_region.name
                    # Define boundaries of current region
                    begin, end = max(curr_region.begin, self.begin), min(curr_region.end, self.end)
                    # Skip if current region does not overlap
                    if end < begin:
                        continue
                    # Otherwise, increase the list for current consensus level
                    regions_list.append(self.__class__(
                        accession=accession, name=name,
                        entry=entry, begin=begin, end=end
                    ))
            # Define regions dictionary
            return regions_dict


def get_consensus(key, msa, summary):
    """ Compute consensus on a multiple sequence alignment

    This method takes as input a multiple sequence alignment with shape (n, m)
    and returns an array with shape (m, ) containing the number of times a
    specific residue has been annotated with a given feature

    Params
    ------
    key : tuple
        Tuple containing feature's evidence, name and consensus
    msa : msa_tools.MSA
        Multiple sequence alignment on which LIP consensus must be retrieved
    summary : Summary
        Summary holding information about sequences and LIP predictions over them

    Returns
    -------
    consensus : np.array
        Array of integers specifying how many times each reisdue has
        been annotated with queried feature key
    """
    # Define alignment shape
    n, m = msa.shape
    # Initialize consensus
    full_consensus = np.zeros((m,), dtype=int)
    # Define set of available accession numbers
    available_acc = set(summary.entries.keys())
    # Loop through each row in the alignment
    for k in range(0, n):
        # Retrieve sequence accession, begin and end
        seq_accession = str(msa.accession[k])
        # Skip if current sequence accession is not in summary
        if seq_accession not in available_acc:
            continue
        # Retrieve sub-sequence begin and end
        seq_begin, seq_end = int(msa.begin[k]), int(msa.end[k])
        # Retrieve binary sequence
        seq_residues = msa.r[k, :]
        # Define aligned positions
        is_aligned = ~msa.is_gap(seq_residues)
        # Get entry in summary for current sequece
        curr_entry = summary.entries.get(seq_accession)
        # Create a feature on current entry
        curr_feature = Summary.Feature(entry=curr_entry, begin=seq_begin, end=seq_end)
        # Retrieve consensus for current feature
        curr_consensus = np.array(curr_feature.get_consensus(key) > 0, dtype=int)
        # Apply consensus on current sequence
        full_consensus[is_aligned] += curr_consensus
    # Return full consensus
    return full_consensus


def get_lip_consensus(*args, **kwargs):
    # Initialize available feature source
    feature_sources = Summary.Feature.feature_sources
    # Define dictionary mapping feature source to consensus
    lip_consensus = dict()
    # Loop through each feature source
    for source_name, source_code in feature_sources.items():
        # Define feature's key for FLIPPER
        feature_id = Summary.Feature.EVIDENCE_DERIVED, Summary.Feature.NAME_LIP, source_code
        # Retrieve consensus for current feature
        lip_consensus.setdefault(source_name, get_consensus(feature_id, *args, **kwargs))
    # Return consensus arrays
    return lip_consensus


def get_disorder_consensus(*args, **kwargs):
    # Initialize available feature source
    feature_sources = Summary.Feature.feature_sources
    # Define dictionary mapping feature source to consensus
    disorder_consensus = dict()
    # Loop through each feature source
    for source_name, source_code in feature_sources.items():
        # Define feature's key for FLIPPER
        feature_id = Summary.Feature.EVIDENCE_CURATED, Summary.Feature.NAME_DISORDER, source_code
        # Retrieve consensus for current feature
        disorder_consensus.setdefault(source_name, get_consensus(feature_id, *args, **kwargs))
    # Return consensus arrays
    return disorder_consensus


# Main
if __name__ == '__main__':

    # Define path to summary file
    path_to_summary = './data/mobidb_result_2021-06-15T15 44 40.082Z.json'
    # Initialize timer
    timer = time.time()
    # Initialize summary as an empty list of entries
    summary = list()
    # Open file
    with open(path_to_summary) as file:
        # Load JSON object
        summary = json.load(file)

    # Print time taken
    print('It took %.02f seconds to load summary file' % (time.time() - timer))
    # Print type of loaded summary
    print('Loaded summary is of type %s' % type(summary))
    # Print space used by loaded summary
    print('Loaded summary has size %.02f MB' % (sys.getsizeof(summary) / 1e06))

    # Get frist entrytuple
    entry = summary[0]
    # Print type of first entry
    print('First entry is of type %s' % type(entry))

    # Loop through each key, value pair in  first entry
    for i, (key, value) in enumerate(entry.items()):
        # Print key and type of value
        print('%d-th item has key %s and value of type %s' % (i, key, type(value)))

    # Print accession
    print('Attribute accession is %s' % entry.get('accession'))
    # Print taxonomy
    print('Attribute taxonomy is %s' % str(entry.get('taxonomy')))
    # Print localization
    print('Attribute localization is %s' % str(entry.get('localization')))

    # Initialize number of entries with full `accession` for accession
    num_full_acc = 0
    # Initialize number of entries with abbreviated `acc` for accession
    num_abbr_acc = 0
    # Define number of entries
    num_entries = len(summary)
    # Loop through each entry
    for entry in summary:
        # Get full `accession` for accession
        num_full_acc += int(entry.get('accession', None) is not None)
        # Get abbreviated `acc` for accession
        num_abbr_acc += int(entry.get('acc', None) is not None)
    # Print number of full accessions
    print('There are %d (%.02f%%) entries with full `accession` attribute' % (
    num_full_acc, num_full_acc / num_entries * 100))
    # Print number of abbreviated accessions
    print('There are %d (%.02f%%) entries with abbreviated `acc` attribute' % (
    num_abbr_acc, num_abbr_acc / num_entries * 100))

    # Retrieve first entry again
    first = summary[0]
    # Load derived LIP predictions
    curated_lip_merge = first.get('derived-lip-merge')
    # Print entry key, value pairs
    print('Values in `curated-lip-merge`')
    for i, (key, value) in enumerate(curated_lip_merge.items()):
        print('  Key %s has value %s' % (key, str(value)))

    # Initialize list of entries which have more than one LIP region
    many_pfam_domains = list()
    # Initialize list of entries which have no LIP region
    zero_pfam_domains = list()
    # Loop through each entry
    for entry in summary:
        # Get Pfam domains
        homology_domain_pfam = entry.get('homology-domain-pfam', dict())
        # Get number of regions
        num_regions = len(homology_domain_pfam.get('regions', list()))
        # Update lists
        if num_regions > 1:
            many_pfam_domains.append(entry)
        elif num_regions < 1:
            zero_pfam_domains.append(entry)
    # Print statistics
    print('There are %d (%.02f%%) entries with more than one Pfam domain' % (
    len(many_pfam_domains), len(many_pfam_domains) / len(summary) * 100))
    print('There are %d (%.02f%%) entries with no Pfam domain' % (
    len(zero_pfam_domains), len(zero_pfam_domains) / len(summary) * 100))

    # Get one entry from those with more than one region
    entry = many_pfam_domains[0]
    # Load derived Pfam mappings
    homology_domain_pfam = entry.get('homology-domain-pfam')
    # Print entry key, value pairs
    print('Values in `homology-domain-pfam`')
    for i, (key, value) in enumerate(homology_domain_pfam.items()):
        print('  Key %s has value %s' % (key, str(value)))

