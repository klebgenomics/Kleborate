"""
Copyright 2025 Mary Maranga (gathonimaranga@gmail.com)
https://github.com/klebgenomics/Kleborate

This file is part of Kleborate. Kleborate is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Kleborate is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Kleborate. If
not, see <https://www.gnu.org/licenses/>.
"""

import os
import json
import re
import pandas as pd
import logging
from ...shared.alignment import align_query_to_ref
LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG) 


# Default argument values 
min_identity=90.0
min_coverage=80.0
percentIdentityOtype=90.0
percentIdentityHtype=95.0
percentCoverageOtype=90.0
percentCoverageHtype=50.0
HIGH_SIMILARITY_THRESHOLD_O=0.00771
MIN_O_IDENTITY_LS=95.0
MIN_O_COVERAGE_LS=48.0


# load database in JSON file
def load_json_file(ref_db):
    """
    Load a JSON file.

    Parameters:
    ref_file  (str): The path to the JSON file.

    Returns:
    dict: The loaded JSON data as a dictionary, or None if an error occurs or the file does not exist.
    """
    if os.path.exists(ref_db):
        try:
            # Open and load the JSON file
            with open(ref_db, 'r') as fp:
                ectyperdb_dict = json.load(fp)
                return ectyperdb_dict
        except Exception as e:
            print(f"An error occurred while loading the JSON file: {e}")
    else:
        print(f"File not found: {ref_db}")
    
    return None


def alignment(ref_file, minimap2_index, assembly, min_identity, min_coverage):
    """
    Aligns the assembly to the reference

    Parameters:
    - assembly: Assembly in FASTA format.
    - ref_file: alleles in FASTA format.
    - minimap2_index: Path to the assembly's minimap2 index
    - min_coverage: Minimum query coverage for alignment.
    - min_identity: Minimum identity percentage for alignment.

    Returns:
    - List of alignment hits in PAF format.
    """

    alignment_hits = align_query_to_ref(
        ref_file,
        assembly,
        ref_index=minimap2_index,
        min_identity=min_identity,
        min_query_coverage=min_coverage
    )
    
    return alignment_hits



def minimap2_output_to_df(alignments):
    """
    Convert the Minimap2 alignment object to a pandas DataFrame.
    
    :param  alignment_hits: List of Minimap2 alignment hits
    :return: DataFrame of the alignment results
    """
    
    # List to store output data
    output_data = []
    
    # Loop over each alignment hit
    for hit in alignments:
        entry = {
            'query_name': hit.query_name,
            'query_length': hit.query_length,
            'query_start': hit.query_start,
            'query_end': hit.query_end,
            'query_cov':hit.query_cov,
            'ref_name': hit.ref_name,
            'ref_length': hit.ref_length,
            'ref_start': hit.ref_start,
            'ref_end': hit.ref_end,
            'percent_identity': hit.percent_identity,
            'alignment_score': hit.alignment_score,
            'ref_seq': hit.ref_seq,
            'strand':hit.strand,
            'ref_cov':hit.ref_cov
        }
        output_data.append(entry)
    
    # Convert the list to a pandas DataFrame
    df = pd.DataFrame(output_data)
    
    return df



def ectyper_dict_to_df(ectyper_dict):
    """
    Convert the ECTyper dictionary into a DataFrame.
    
    :param ectyper_dict: ECTyper database as a dictionary from a JSON file of all known O and H alleles.
    :return: DataFrame containing the information from the JSON file.
    """

    # Initialize an empty list to store the records that will be converted to a DataFrame.
    temp_list = []

    # Define the antigen types to look for in the dictionary: O and H antigens.
    antigens = ["O", "H"]

    # Loop through each antigen type (O and H).
    for antigen in antigens:
        # Retrieve the alleles for the current antigen type.
        alleles = ectyper_dict[antigen]

        # Iterate over each allele entry in the dictionary.
        for name, allele in alleles.items():
            # Create a dictionary for each allele with its relevant information.
            new_entry = {
                'type': antigen,  # Antigen type (O or H)
                'antigen': allele.get('allele'),  # Antigen allele
                'name': name,  # Name of the allele
                'gene': allele.get('gene'),  # Gene name
                'desc': allele.get('desc'),  # Description of the allele
                'sharedallele': allele.get('isAlleleShared')  # Whether the allele is shared
            }
            # Append the new entry to the list.
            temp_list.append(new_entry)
    
    # Convert the list of dictionaries to a pandas DataFrame.
    df = pd.DataFrame(temp_list)

    # Return the DataFrame for further use.
    return df



def setAlleleMeta(keys,metadf):
    metaalleledict = {}
    for key in keys:
        metaalleledict[key] = metadf[key]
    return metaalleledict



def get_prediction(per_genome_df, args):
    
    """
    Make serotype prediction for a single genome based on the minimap2 output.

    :param per_genome_df: The minimap2 reference allele hits results for a given genome.
    :param args: Commandline arguments.
    :return: Serotype dictionary with O and H antigen information.
    """
    # The per_genome_df DataFrame is sorted in descending order based on the 'score' column.
    per_genome_df = per_genome_df.sort_values(by=['alignment_score'], ascending=False, kind='mergesort')
    # Defining the dictionary keys that will hold serogroup and genescores.
    serotype_dictkeys = {"serogroup", "genescores"}

    # Initializing the serotype dictionary for O and H antigens with default values as '-'.
    serotype = {
        'O': dict.fromkeys(serotype_dictkeys, "-"),
        'H': dict.fromkeys(serotype_dictkeys, "-")
    }

    # Initializing genescores and alleles for both O and H antigens.
    serotype["O"]["genescores"] = {}
    serotype["O"]["alleles"] = {}
    serotype["H"]["genescores"] = {}
    serotype["H"]["alleles"] = {}

    # Go for the highest match, if both genes exist over the thresholds
    minimapresultsdict={}; minimapresultsdict["O"]={}; minimapresultsdict["H"]={}

    # Iterating over each row of the DataFrame (per_genome_df).
    for row in per_genome_df.itertuples():
        # H is already set, skip
        # get the 'O' or 'H' from the antigen column
        ant = row.antigen[:1]
        if row.query_name not in minimapresultsdict[ant]:  #row.qseqid = allele database key (e.g. H29-2-fliC-origin)
            minimapresultsdict[ant][row.query_name] = {}
            minimapresultsdict[ant][row.query_name]["gene"] = row.gene
            minimapresultsdict[ant][row.query_name]["antigen"] = row.antigen
            minimapresultsdict[ant][row.query_name]["score"] = float(row.alignment_score)
            minimapresultsdict[ant][row.query_name]["identity"] = float(row.percent_identity)
            minimapresultsdict[ant][row.query_name]["coverage"] = float(row.query_cov)# used query_cov instead of qcovhsp
            minimapresultsdict[ant][row.query_name]["coverage"] = float(row.ref_cov)
            minimapresultsdict[ant][row.query_name]["contigname"] = row.ref_name
            minimapresultsdict[ant][row.query_name]["startpos"] = int(row.ref_start)
            minimapresultsdict[ant][row.query_name]["endpos"] = int(row.ref_end)
            minimapresultsdict[ant][row.query_name]["length"] = int(row.ref_length)
            minimapresultsdict[ant][row.query_name]["shared"] = bool(row.sharedallele)
    allelefieldnames = ["identity", "coverage", "contigname", "length", "startpos", "endpos", "gene"]

    ant = "H"
    if minimapresultsdict[ant].keys():
        topHallele = sorted([(dballele,minimapresultsdict[ant][dballele]["score"]) for dballele in minimapresultsdict["H"].keys()],
                            key=lambda x:x[1], reverse=True)[0][0]
        serotype[ant]["serogroup"] = minimapresultsdict[ant][topHallele]["antigen"]
        gene = minimapresultsdict[ant][topHallele]["gene"]
        score = minimapresultsdict[ant][topHallele]["score"]
        serotype[ant]["genescores"] = {gene:score}
        serotype[ant]["alleles"][topHallele] = setAlleleMeta(allelefieldnames,minimapresultsdict[ant][topHallele])
    else:
        LOG.warning("No H antigen alleles were identified in this sample")

    sortedOalleles = [tuple[0] for tuple in sorted([(dballele,
                                             minimapresultsdict["O"][dballele]["score"])
                                             for dballele in minimapresultsdict["O"].keys()],
                                             key=lambda x: x[1], reverse=True)]

    otype={}
    for allele in sortedOalleles:
        oantigen = minimapresultsdict["O"][allele]["antigen"]
        if minimapresultsdict["O"][allele]["antigen"] not in otype.keys():
            otype[oantigen] = {"genescores":{}, "alleles":[], "allele2gene":{}}
        if minimapresultsdict["O"][allele]["gene"] not in otype[oantigen]["genescores"].keys():
            gene = minimapresultsdict["O"][allele]["gene"]
            otype[oantigen]["genescores"][gene] = minimapresultsdict["O"][allele]["score"]
            otype[oantigen]["alleles"].append(allele)
            otype[oantigen]["allele2gene"][allele] = gene

    # rank O-type serovars based on the sum of scores from BOTH alleles (wzx/wzy or wzt/wzm)
    # calculate score  for O antigen
    rank_Otype_dict = {}
    for oantigen in otype.keys():
        rank_Otype_dict[oantigen]={"numalleles": 0, "scores":[], "sumscore": 0}
        rank_Otype_dict[oantigen]["numalleles"] = len(otype[oantigen]["alleles"])
        if rank_Otype_dict[oantigen]["numalleles"] == 4:
            LOG.warning("O-antigen {} has 4 alleles instead of expected 2."
                        "Might be due extra genes pairs (wzx/wzy,wzt/wzm) mapping to the same antigen or multiple gene copies".format(oantigen))
            wzx_wzy_score_sum = otype[oantigen]["genescores"]['wzx'] + otype[oantigen]["genescores"]['wzy']
            wzm_wzt_score_sum = otype[oantigen]["genescores"]['wzm'] + otype[oantigen]["genescores"]['wzt']
            if wzx_wzy_score_sum > wzm_wzt_score_sum:
                rank_Otype_dict[oantigen]["scores"] = [otype[oantigen]["genescores"][gene] for gene in ['wzx','wzy']]
            else:
                rank_Otype_dict[oantigen]["scores"] = [otype[oantigen]["genescores"][gene] for gene in ['wzm', 'wzt']]

        elif rank_Otype_dict[oantigen]["numalleles"] == 3:
            if 'wzx' in otype[oantigen]["genescores"].keys() and 'wzy' in otype[oantigen]["genescores"].keys():
                rank_Otype_dict[oantigen]["scores"] = [otype[oantigen]["genescores"][gene] for gene in ['wzx', 'wzy']]
            elif 'wzm' in otype[oantigen]["genescores"].keys() and 'wzt' in otype[oantigen]["genescores"].keys():
                rank_Otype_dict[oantigen]["scores"] = [otype[oantigen]["genescores"][gene] for gene in ['wzm', 'wzt']]
        else:
            rank_Otype_dict[oantigen]["scores"] = otype[oantigen]["genescores"].values()
        rank_Otype_dict[oantigen]["sumscore"] = sum(rank_Otype_dict[oantigen]["scores"])

    scorestupleslist = [(otypename,rank_Otype_dict[otypename]["sumscore"]) for otypename in rank_Otype_dict.keys()]
    scorestupleslist = sorted(scorestupleslist, key=lambda x: x[1], reverse=True) #[('O102', 1.73)]
    best_order_list = [item[0] for item in scorestupleslist] #['O102']


    LOG.debug("Otype dict:{}".format(otype))
    LOG.debug("Serotype dict:{}".format(serotype))
    LOG.debug("Best order alleles-scores list of tuples:{}".format(scorestupleslist))
    LOG.debug("Best order alleles list:{}".format(best_order_list))

    # having gone through all the hits over the threshold, make the call
    # go through the O-antigens in order, making the call on the first that have
    # a matching pair

    ant="O"; selectedOantigen = ""
    for oantigen in best_order_list:
        # if wzm / wzy or wzx / wzy, call the match
        if 'wzx' in otype[oantigen]["genescores"].keys() and 'wzy' in otype[oantigen]["genescores"].keys():
            serotype[ant]['serogroup'] = oantigen
            serotype[ant]['genescores'] = {'wzx':otype[oantigen]["genescores"]["wzx"],
                                           'wzy':otype[oantigen]["genescores"]["wzy"]}
            selectedOantigen = oantigen
            break
        elif 'wzm' in otype[oantigen]["genescores"].keys() and 'wzt' in otype[oantigen]["genescores"].keys():
            serotype[ant]['serogroup'] = oantigen
            serotype[ant]['genescores'] = {'wzm': otype[oantigen]["genescores"]["wzm"],
                                           'wzt': otype[oantigen]["genescores"]["wzt"]}
            selectedOantigen = oantigen
            break

        # FIX: O-antigen typing might fail due to poor sequencing or inability to assemble one of the wzx/wzy/wzm/wzt loci
        # if only one of the signatures is found, still produce output but warn user on false positives
        elif 'wzx' in otype[oantigen]["genescores"].keys() or 'wzy' in otype[oantigen]["genescores"].keys():
            serotype['O']['serogroup']  = oantigen

            if 'wzx' in otype[oantigen]["genescores"].keys():
                serotype[ant]['genescores'] = {'wzx': otype[oantigen]["genescores"]["wzx"]}
            elif 'wzy' in otype[oantigen]["genescores"].keys():
                serotype[ant]['genescores'] = {'wzy': otype[oantigen]["genescores"]["wzy"]}
            selectedOantigen = oantigen
            break

        elif 'wzm' in otype[oantigen]["genescores"].keys() or 'wzt' in otype[oantigen]["genescores"].keys():
            serotype['O']['serogroup'] = oantigen
            if 'wzm' in otype[oantigen]["genescores"].keys():
                serotype[ant]['genescores'] = {'wzm': otype[oantigen]["genescores"]["wzm"]}
            elif 'wzt' in otype[oantigen]["genescores"].keys():
                serotype[ant]['genescores'] = {'wzt': otype[oantigen]["genescores"]["wzt"]}
            selectedOantigen = oantigen
            break

    identicalscorestupleslist = [(orow, ocol, abs(i - j)) for ocol, i in scorestupleslist for orow, j in scorestupleslist
                            if i - j == 0 and (orow == selectedOantigen or ocol == selectedOantigen) and
                            (orow != ocol)]

    #find highly similar antigens at 99.9% indentity which is only 1 nt apart
    highsimilarity_oantigens = [orow for ocol, i in scorestupleslist for orow, j in
                                scorestupleslist if i - j <= HIGH_SIMILARITY_THRESHOLD_O
                                and orow != selectedOantigen
                                and (orow == selectedOantigen or ocol == selectedOantigen)
                                and (orow != ocol)]

    if highsimilarity_oantigens:
        mixedoantigen = [selectedOantigen] + highsimilarity_oantigens
        serotype['O']['serogroup'] = "/".join(mixedoantigen)
        LOG.info("Highly similar O-antigen candidates were found for {}".format(mixedoantigen))
    elif selectedOantigen != "":
        serotype['O']['serogroup'] = selectedOantigen

    # append to existing top O-antigen and generate mixed antigen call
    # alleles with identical score
    if identicalscorestupleslist:
        identical_oantigens = [i for i, j, s in identicalscorestupleslist if i != selectedOantigen]

        for oantigen in identical_oantigens:
            for allele in otype[oantigen]["alleles"]:
                serotype[ant]["alleles"][allele] = setAlleleMeta(allelefieldnames, minimapresultsdict[ant][allele])


    # add info for selected alleles
    if selectedOantigen:
        for allele in otype[selectedOantigen]["alleles"]:
            serotype[ant]["alleles"][allele] = setAlleleMeta(allelefieldnames, minimapresultsdict[ant][allele])

    return serotype


def predict_serotype(alignment_hits, ectyper_dict, args, previous_results):
    """
    Predict the serotype given the minimap output of the markers against the genomes

    :param alignment_hits: minimap results of O and H-type allele search against the genomes of interest
    :param ectyper_dict: ectyper database in dict format from JSON file of known alleles and their O and H mappings
    :param args: Commandline arguments
    :return: The CSV formatted predictions file
    """

    LOG.info("Predicting serotype from minimap output")
    output_df = minimap2_output_to_df(alignment_hits)

    ectyper_df = ectyper_dict_to_df(ectyper_dict) 

    # Merge output_df and ectyper_df
    output_df = output_df.merge(ectyper_df, left_on='query_name', right_on='name', how='left')
    predictions_dict = {}

    # Get potential O antigens
    OantigensPotential = set(output_df.query('type == "O" & percent_identity > ' + str(MIN_O_IDENTITY_LS) +
                                             ' & query_cov > ' + str(MIN_O_COVERAGE_LS))["antigen"])

    # Filter the output_df based on different conditions
    if len(OantigensPotential) == 1:
        output_df = output_df.query(
            '(type == "O" & percent_identity >= ' + str(MIN_O_IDENTITY_LS) + ' & query_cov  >= ' +
            str(MIN_O_COVERAGE_LS) + ') | '
            '(type == "H" & percent_identity >= ' + str(percentIdentityHtype) + ' & query_cov  >= ' + str(percentCoverageHtype) + ' )'
        )
    else:
        output_df = output_df.query(
            '(type == "O" & percent_identity >= ' + str(percentIdentityOtype) + ' & query_cov  >= ' + str(percentCoverageOtype) + ') | '
            '(type == "H" & percent_identity >= ' + str(percentIdentityHtype) + ' & query_cov  >= ' + str(percentCoverageHtype) + ' )'
        )

    # Assign genome name from previous_results
    output_df = output_df.assign(genome_name=previous_results['strain'])

    # Make prediction for each genome based on minimap2 output
    for genome_name, per_genome_df in output_df.groupby('genome_name'):
        predictions_dict[genome_name] = get_prediction(per_genome_df, args)

    return predictions_dict, output_df



def mean(numbers):
    return sum(numbers)/len(numbers)


def getPredictionNumAlleles(sample, final_results_dict):
    numalleles = 0
    if "O" in final_results_dict[sample]:
        numalleles += len(final_results_dict[sample]["O"]["genescores"].keys()) #some O-antigens like O8 have both wzx/wzy and wzm/wzt
    if "H" in final_results_dict[sample]:
        numalleles += len(final_results_dict[sample]["H"]["genescores"].keys())
    return numalleles


def getQuality_control_results(sample, final_results_dict, ectyperdb_dict):
    """
    Determined approximate quality of the prediction based on the allele scores
    :param sample: sample/genome name)
    :param final_results_dict: dictionary with final output results (e.g. serovar, sequences, conf. scores)
    :return: Quality dictionary
    """
    
    # Initialize error key if it does not exist
    if "error" not in final_results_dict[sample]:
        final_results_dict[sample]["error"] = ""

    if 'O' in final_results_dict[sample]:
        Otype = final_results_dict[sample]['O']['serogroup']
        Otypealleles = final_results_dict[sample]['O']['alleles'].keys()
    else:
        Otype = "-"
        Otypealleles = []

    if 'H' in final_results_dict[sample]:
        Htype = final_results_dict[sample]['H']['serogroup']
        Htypealleles = final_results_dict[sample]['H']['alleles'].keys()
    else:
        Htype = "-"
        Htypealleles = []

    if not 'O' and 'H' in final_results_dict[sample]:
        return "-"

    # QC2: Check on quality of O and H antigen prediction if at all taking into account resolution power of each antigen
    if Otype == "-" and Htype == "-":
        final_results_dict[sample]["error"] += "Failed to type both O and H antigens. Consider lowering thresholds or traditional serotyping. "
        return "FAIL (-:- TYPING)"
    elif Otype == "-" and Htype != "-":
        final_results_dict[sample]["error"] += "Failed to type O antigen. Consider lowering thresholds or traditional serotyping. "
        return "WARNING (-:H TYPING)"
    elif Otype != "-" and Htype == "-":
        final_results_dict[sample]["error"] += "Failed to type H antigen. Consider lowering thresholds or traditional serotyping. "
        return "WARNING (O:- TYPING)"
    elif len(Otype.split("/")) >= 2:
        final_results_dict[sample]["error"] += "Mixed O-antigen call reported due to very high degree of antigen groups identity in excess of {}%. Consider traditional serotyping as in-silico predictions might not be accurate.".format((1 - HIGH_SIMILARITY_THRESHOLD_O) * 100)
        return "WARNING MIXED O-TYPE"
    else:
        checkpredscoresboolO = [final_results_dict[sample]["O"]["alleles"][allele]["identity"] *
                                final_results_dict[sample]["O"]["alleles"][allele]["coverage"] / 1e4 >=
                                ectyperdb_dict["O"][allele]["MinPident"] * ectyperdb_dict["O"][allele]["MinPcov"] / 1e4
                                for allele in Otypealleles]

        checkpredscoresboolH = [final_results_dict[sample]["H"]["alleles"][allele]["identity"] *
                                final_results_dict[sample]["H"]["alleles"][allele]["coverage"] / 1e4 >=
                                ectyperdb_dict["H"][allele]["MinPident"] * ectyperdb_dict["H"][allele]["MinPcov"] / 1e4
                                for allele in Htypealleles]

        if all(checkpredscoresboolO + checkpredscoresboolH) and all([item != [] for item in [checkpredscoresboolO, checkpredscoresboolH]]):
            QCflag = "PASS (REPORTABLE)"
        elif all(checkpredscoresboolO):
            QCflag = "WARNING (H NON-REPORT)"
            final_results_dict[sample]["error"] += "H-antigen has %identity or %coverage below reportable threshold. "
        elif all(checkpredscoresboolH):
            QCflag = "WARNING (O NON-REPORT)"
            final_results_dict[sample]["error"] += "O-antigen has %identity or %coverage below reportable threshold. "
        else:
            QCflag = "WARNING (O and H NON-REPORT)"
            final_results_dict[sample]["error"] += "O and H-antigens have %identity or %coverage below reportable threshold. "

    return QCflag


def process_sample_data(sample, final_dict, ectyperdb_dict):
    # Add Quality Control Flag (QC)
    final_dict[sample]["QC"] = getQuality_control_results(sample, final_dict, ectyperdb_dict)

    # Extract serogroup information for O and H antigens
    Otype = final_dict[sample]["O"]["serogroup"] if "O" in final_dict[sample] else "-"
    Htype = final_dict[sample]["H"]["serogroup"] if "H" in final_dict[sample] else "-"
    
    # Combine Otype and Htype to form Serotype
    Serotype = f"{Otype}:{Htype}" if Otype != "-" and Htype != "-" else "-"

    # Initialize strings and lists for gene scores and alleles
    genescoresstr = ""
    geneordertuplelist = []
    genelistordersimple = []
    alleles_list = []

    # Process O genescores
    num_alleles = 0
    if Otype != "-":
        for Ogene, score in sorted(final_dict[sample]["O"]["genescores"].items()):
            genescoresstr += "{}:{:.3f};".format(Ogene, score)  # No scientific notation, keeping decimal places
            geneordertuplelist.append(("O", Ogene, score))
            genelistordersimple.append(Ogene)
            num_alleles += 1

        # Add alleles for O genes
        if "alleles" in final_dict[sample]["O"]:
            alleles_list.extend(final_dict[sample]["O"]["alleles"].keys())

    # Process H genescores
    if Htype != "-":
        for Hgene, score in final_dict[sample]["H"]["genescores"].items():
            genescoresstr += "{}:{:.3f};".format(Hgene, score)  # No scientific notation, keeping decimal places
            geneordertuplelist.append(("H", Hgene, score))
            genelistordersimple.append(Hgene)

        # Add alleles for H genes
        if "alleles" in final_dict[sample]["H"]:
            alleles_list.extend(final_dict[sample]["H"]["alleles"].keys())

    # Create the O allele count information (renamed from output_line to O_allele_No)
    O_allele_No = []
    if num_alleles == 0:
        O_allele_No.append("-")
    else:
        O_allele_No.append("Based on {} allele(s)".format(num_alleles))

    # Append gene scores string or "-" if empty
    if genescoresstr == "":
        O_allele_No.append("-")
    else:
        O_allele_No.append(genescoresstr)

    # Prepare individual report fields for allele information
    report_data = {}
    for ant, gene, score in geneordertuplelist:
        # Iterate through the alleles to find the correct one
        if "alleles" in final_dict[sample][ant]:
            found_allele = False
            for allele_key, allele_info in final_dict[sample][ant]["alleles"].items():
                if allele_info["gene"] == gene:  # Match based on gene name
                    report_data[gene] = {
                        "identity": allele_info.get("identity", "-"),
                        "coverage": allele_info.get("coverage", "-"),
                        "contigname": allele_info.get("contigname", "-"),
                        "coordinates": f"{allele_info.get('startpos', '-')}-{allele_info.get('endpos', '-')}" if "startpos" in allele_info and "endpos" in allele_info else "-",
                        "length": allele_info.get("length", "-")
                    }
                    found_allele = True
                    break  # Stop searching once the correct allele is found
            # If no allele is found for the gene, set empty values
            if not found_allele:
                report_data[gene] = {
                    "identity": "-",
                    "coverage": "-",
                    "contigname": "-",
                    "coordinates": "-",
                    "length": "-"
                }

    # Compile all information into a single results dictionary
    results = {
        "Serotype": Serotype,
        "Otype": Otype,
        "Htype": Htype,
        "gene_scores": genescoresstr,
        "O_allele_No": O_allele_No,
        "alleles": ";".join(alleles_list),
        "QC": final_dict[sample]["QC"],
        **report_data
    }

    return results

