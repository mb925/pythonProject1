import sys
import config as cfg
import dam_script as dm
import json
import pandas as pd


def main():
    # res = dm.Summary.load_ndjson(cfg.data['mobi'] + '/mobidb.json')
    #
    # print(res)

    # df_repeatsdb_units = pd.DataFrame(columns=["PDB", "type", "start", "end",  "review", "annotator", "origin"])
    # df_repeatsdb_units.to_csv(cfg.data['repeatsdb'] + '/repeatsdb_units_ins.tsv', index=False, sep='\t')
    #
    #
    # # entries json to csv
    # with open(cfg.data['repeatsdb'] + '/repeatsdb_units_ins.tsv', 'a') as f:
    #     for line in open(cfg.data['repeatsdb'] + '/entries.mjson', 'r'):
    #             line = json.loads(line)
    #             rev = '0'
    #             if line['reviewed']:
    #               rev = '1'
    #
    #             f.write(line['pdb_id'] + line['pdb_chain'] + "\t" +
    #                     line['type'] + "\t" + line['start'] + "\t" + line['end'] + "\t" +
    #                     rev + "\t" +  line['annotator'] + "\t" + line['origin'] + "\n")

    # filter repeatsdb
    df_repeatsdb = pd.read_csv(cfg.data['repeatsdb'] + '/RepeatsDB-table.tsv', sep='\t')
    # filter repeatsdb table (Homo Sapiens)
    df_repeatsdb = df_repeatsdb[df_repeatsdb['Title'].str.contains("Homo")]
    df_repeatsdb_uniprot = df_repeatsdb.set_index('UniProt')['PDB chain IDs'].str.split(',', expand=True).stack().reset_index('UniProt')
    df_repeatsdb_uniprot.rename(columns={0: 'PDB'}, inplace=True)
    # adding uniprot info
    df_repeatsdb_units = pd.read_csv(cfg.data['repeatsdb'] + '/repeatsdb_units_ins.tsv', sep='\t')
    print(len(df_repeatsdb_units))
    df_repeatsdb_units_uniprot = df_repeatsdb_units.merge(df_repeatsdb_uniprot, how='inner', on='PDB')
    df_repeatsdb_units_uniprot.to_csv(cfg.data['repeatsdb'] + '/df_repeatsdb_units_human_uniprot.tsv', index=False, sep='\t')
    print(df_repeatsdb_uniprot)

if __name__ == '__main__':
    sys.exit(main())
