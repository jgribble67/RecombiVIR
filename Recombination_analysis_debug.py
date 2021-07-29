#!/bin/python3
##Last modified 07/14/21 by JGB
# import argparse
import pandas as pd
import os
import fnmatch
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import math

# parser = argparse.ArgumentParser()
# parser.add_argument("Sample_list", help="A text file with each sample base name on a new line.")
# parser.add_argument("Virus", help="Virus name. Options are MHV, MERS, SARS2.")
# parser.add_argument("Working_Dir", help="Absolute or relative path of directory with data.")
# parser.add_argument("Experiment_Name", help="Experiment name for naming output reports.")
# parser.add_argument("--verison", help="Version of ViReMa utilized. Default is 0.21")
# parser.add_argument("--Output_Dir", help="Absolute or relative path of directory for output folders and files. Default working directory.")
# parser.add_argument("--Shannon_Entropy", help="Path to folder with Virus_Recombination_Results.txt files for Shannon Entropy")
# parser.add_argument("--Virus_Accession", help="NCBI virus accession number")
# parser.add_argument("--Min_Coverage", help="Minimum counts to include in calculation of Shannon Entropy")
# args = parser.parse_args()

#Make a report dataframe with sample column loaded
version = 0.20
if (version >= 0.21):
    report = pd.DataFrame(columns=['sample',
                                   'unique_junctions',
                               'recombined_nts',
                               "total_nts",
                               "total_cutting_f_nts",
                               "total_cutting_r_nts",
                               "total_cutting_site_nts",
                               "cutting_f_jfreq",
                               "cutting_r_jfreq",
                               "cutting_jfreq",
                               "jfreq"])
if (version < 0.21):
    report = pd.DataFrame(columns=['sample',
                                   'unique_junctions',
                                   'recombined_nts',
                                   "total_nts",
                                   "jfreq"])
# sample_list = [line.rstrip('\n') for line in open(str(args.Sample_list))]
sample_list = [line.rstrip('\n') for line in open("/home/denison-thelio/Current_projects/RNAseq/MA_MHV_NHC_virion/samples.txt")]
report['sample'] = sample_list
#Set other variables
virus = "MHV"
wd = "/home/denison-thelio/Current_projects/RNAseq/MA_MHV_NHC_virion/"
od = wd
exp = "MA_MHV_NHC_virion"

# virus = str(args.Virus)
# wd = str(args.Working_Dir)
# if args.Output_Dir:
#     od = str(args.Output_Dir)
# else:
#     od = wd
# exp = str(args.Experiment_Name)
# if args.version:
#     version = float(args.version)
# else:
#     version = 0.21

#Make target folders for output files
if not os.path.exists(od + '/Junction_Files'):
    os.makedirs(od + '/Junction_Files')
save_dir_file = od + 'Junction_Files/'
if not os.path.exists(od + '/Junction_Plots'):
    os.makedirs(od + '/Junction_Plots')
save_dir_plot = od + "Junction_Plots/"

##Shannon Entropy script originally authored by Andrew Routh.
Shannon_Entropy = True
if Shannon_Entropy == True:
    Min_Coverage = 0
    # if args.Min_Coverage:
    #     Min_Coverage = int(args.Min_Coverage)
    # else:
    #     Min_Coverage = 0
    # se_dir = str(args.Shannon_Entropy)
    se_dir = od + "Shannon_Entropy/"
    # Virus_Accession = str(args.Virus_Accession)
    Virus_Accession = "AY910861.1"
    se_output_normalized = pd.DataFrame(columns=['sample', Virus_Accession + "_to_" + Virus_Accession, Virus_Accession + "_RevStrand_to_" + Virus_Accession, Virus_Accession + "_RevStrand_to_" + Virus_Accession + "_RevStrand", Virus_Accession + "_to_" + Virus_Accession + "_RevStrand"])
    # se_output_normalized = pd.DataFrame(columns=['sample'])
    se_output_normalized['sample'] = sample_list
    # se_output = pd.DataFrame(columns=['sample'])
    # se_output['sample'] = sample_list
    for file in os.listdir(se_dir):
        if fnmatch.fnmatch(file, "*_Virus_Recombination_Results.txt"):
            sample_name = str(file.split("_")[0])
            Dicts = {}
            with open(se_dir + file, 'r') as file1:
                Data = file1.readline()
                while Data:
                    Name = Data[13:-1]
                    Dicts[Name] = file1.readline().split("\t")[:-1]
                    Data = file1.readline()
                    Data = file1.readline()
            DictKeys = {}
            n = 1
            for Gene in Dicts:
                Data = Dicts[Gene]
                if Virus_Accession in Gene:
                    coverage_file = pd.read_csv(wd + sample_name + "_virema/" + sample_name + "_virema_coverage.txt", sep = "\t", header=0)
                    Virus_Coverage = np.mean(coverage_file['Coverage'])
                    Total_Reads = Virus_Coverage
                else:
                    Total_Reads = 0
                    print("Running Shannon Entropy Calculation for " + sample_name + ". Unknown genome and not normalizing to coverage.")
                Sums = []
                Rec_Total = 0
                for i in Data:
                    data = i.split("_")
                    Freq = int(data[-1])
                    Rec_Total += Freq
                    Sums.append(Freq)
                Entropy = 0
                # for i in Sums:
                #     Fraction = i / float(Rec_Total)
                #     Entropy -= math.log(Fraction, 2) * Fraction
                #     se_output.loc[se_output["sample"] == sample_name, [str(Gene)]] = Entropy
                #     se_output.to_csv(od + sample_name + "_shannon_entropy.txt", sep="\t", index=False)
                # Entropy = 0
                for i in Sums:
                    Fraction = i / float(Rec_Total + Total_Reads)
                    Entropy -= math.log(Fraction, 2) * Fraction
                Fraction = Total_Reads / float(Rec_Total + Total_Reads)
                Entropy -= math.log(Fraction, 2) * Fraction
                se_output_normalized.loc[se_output_normalized["sample"] == sample_name, [str(Gene)]] = Entropy
se_output_normalized['mean_entropy'] = se_output_normalized.mean(axis=1)
se_output_normalized.to_csv(od + exp + "_shannon_entropy_normalized.txt", sep="\t", index=False)
#Isolate forward junctions and make junction plots.
bed_dir = wd + "BED_Files/"
for file in os.listdir(bed_dir):
    if fnmatch.fnmatch(file, "*_Virus_Recombination_Results.bed"):
        sample_name = str(file.split("_")[0])
        if (version >= 0.21):
            bed = pd.read_csv(bed_dir + file, sep="\t", header=0, index_col=False, names=['genome', 'start', 'stop', 'type', 'depth', 'strand', 'rgb1', 'rgb2', 'start_seq', 'stop_seq'])
        if (version < 0.21):
            bed = pd.read_csv(bed_dir + file, sep="\t", header=0, index_col=False,
                              names=['genome', 'start', 'stop', 'type', 'depth', 'strand', 'start1', 'stop1'])
            bed = bed.drop(['start1', 'stop1'], axis=1)
        unique_junctions = len(bed.index)
        recombined_nts = bed['depth'].sum()
        bed = bed.sort_values(by=['depth'], ascending=True)
        total = bed['depth'].sum()
        bed['frequency'] = bed['depth'] / total
        bed['logfreq'] = np.log10(bed['frequency'])
        bed = bed.reset_index(drop=True)
        bed_forward = bed.loc[bed['start'] < bed['stop']]
        bed_forward = bed_forward.reset_index(drop=True)
        report.loc[report['sample'] == str(sample_name), ['unique_junctions']] = unique_junctions
        report.loc[report['sample'] == str(sample_name), ['recombined_nts']] = recombined_nts
        bed.to_csv(save_dir_file + sample_name + '_junctions.txt', sep='\t', index=False)
        bed_forward.to_csv(save_dir_file + sample_name + '_forward_junctions.txt', sep='\t', index=False)
        if (virus == 'MHV'):
            sns.set_style("ticks")
            fig = plt.figure(figsize=(4, 4))
            plt.ioff()
            plt.scatter(bed_forward.stop, bed_forward.start, c=bed_forward.logfreq, cmap='gist_rainbow', alpha=1, vmin=0, vmax=-6, s=15)
            plt.xlim([-1500, 33500])
            plt.ylim([-1500, 33500])
            plt.xticks(fontsize=10)
            plt.yticks(fontsize=10)
            plt.xlabel("3' Positon", fontsize=14)
            plt.ylabel("5' Position", fontsize=14)
            cax = fig.add_axes([0.15, 0.95, 0.70, 0.02])
            cbar = plt.colorbar(orientation="horizontal", cax=cax)
            cbar.ax.tick_params(labelsize=10)
            cbar.ax.set_title("log10(Frequency)", fontsize=12)
            plt.savefig(save_dir_plot + sample_name + "_junctionplot.png", dpi=600, bbox_inches='tight')
            plt.close('all')
        if (virus == 'MERS' or virus == 'SARS2'):
            sns.set_style("ticks")
            fig = plt.figure(figsize=(4, 4))
            plt.ioff()
            plt.scatter(bed_forward.stop, bed_forward.start, c=bed_forward.logfreq, cmap='gist_rainbow', alpha=1,
                        vmin=0, vmax=-6, s=15)
            plt.xlim([-500, 31500])
            plt.ylim([-500, 31500])
            plt.xticks(fontsize=10)
            plt.yticks(fontsize=10)
            plt.xlabel("3' Positon", fontsize=14)
            plt.ylabel("5' Position", fontsize=14)
            cax = fig.add_axes([0.15, 0.95, 0.70, 0.02])
            cbar = plt.colorbar(orientation="horizontal", cax=cax)
            cbar.ax.tick_params(labelsize=10)
            cbar.ax.set_title("log10(Frequency)", fontsize=12)
            plt.savefig(save_dir_plot + sample_name + "_junctionplot.png", dpi=600, bbox_inches='tight')
            plt.close('all')
    if (version >= 0.21):
        if fnmatch.fnmatch(file, "*_Virus_cuttingsites.f.bedgraph"):
            sample_name = str(file.split("_")[0])
            depth_f = pd.read_csv(bed_dir + file, sep="\t", header = 0, index_col=False, names=['genome', 'position', 'position1', 'coverage'])
            f_nts = sum(depth_f['coverage'])
            report.loc[report['sample'] == sample_name, ['total_cutting_f_nts']] = f_nts
        if fnmatch.fnmatch(file, "*_Virus_cuttingsites.r.bedgraph"):
            sample_name = str(file.split("_")[0])
            depth_r = pd.read_csv(bed_dir + file, sep="\t", header=0, index_col=False, names=['genome', 'position', 'position1', 'coverage'])
            r_nts = sum(depth_r['coverage'])
            report.loc[report['sample'] == sample_name, ['total_cutting_r_nts']] = r_nts
for file in os.listdir(wd):
    if fnmatch.fnmatch(file, "*_coverage.txt"):
        sample_name = file.split("_")[0]
        depth = pd.read_csv(wd + file, sep="\t", header = 0)
        total_depth = sum(depth['Coverage'])
        report.loc[report['sample'] == sample_name, ['total_nts']] = total_depth
if (version >= 0.21):
    report['total_cutting_site_nts'] = report['total_cutting_f_nts'] + report['total_cutting_r_nts']
    report['cutting_f_jfreq'] = (report['recombined_nts'] / report['total_cutting_f_nts']) * 1000000
    report['cutting_r_jfreq'] = (report['recombined_nts'] / report['total_cutting_r_nts']) * 1000000
    report['cutting_jfreq'] = (report['recombined_nts'] / report['total_cutting_site_nts']) * 1000000
report['jfreq'] = (report['recombined_nts'] / report['total_nts']) * 1000000
report.to_csv(od + exp + "_ViReMa_report.txt", sep="\t", index=False)

##sgmRNA filtering and quantification
if not os.path.exists(od + '/sgmRNAs_DVGs'):
    os.makedirs(od + '/sgmRNAs_DVGs')
save_dir_sgmRNAs = od + 'sgmRNAs_DVGs/'
sgmRNA_report = pd.DataFrame(columns=['sample',
                                      'total_nts',
                                      'total_sgmRNA_depth',
                                      'major_sgmRNA_depth',
                                      'minor_sgmRNA_depth',
                                      'DVG_depth',
                                      'total_junctions',
                                      'percent_DVGs',
                                      'percent_major_sgmRNA',
                                      'percent_minor_sgmRNA',
                                      'DVG_jfreq',
                                      'major_sgmRNA_jfreq',
                                      'minor_sgmRNA_jfreq'])
sgmRNA_report['sample'] = sample_list
for file in os.listdir(wd):
    if fnmatch.fnmatch(file, "*_coverage.txt"):
        sample_name = file.split("_")[0]
        depth = pd.read_csv(wd + file, sep="\t", header = 0)
        total_depth = sum(depth['Coverage'])
        sgmRNA_report.loc[sgmRNA_report['sample'] == sample_name, ['total_nts']] = total_depth
for file in os.listdir(wd + "Junction_Files/"):
    if fnmatch.fnmatch(file, "*_forward_junctions.txt"):
        sample_name = file.split("_")[0]
        forward_junctions = pd.read_csv(wd + "Junction_Files/" + file, sep="\t", header=0)
        if (virus == "MHV"):
            forward_junctions['start_type'] = forward_junctions['start'].apply(
                lambda x: "TRSL" if ((x >= 32) & (x <= 102)) else "DVG")
            forward_junctions['stop_type'] = forward_junctions['stop'].apply(
                lambda x: "sgmRNA2" if ((x >= 21714) & (x <= 21784)) else (
                    "sgmRNA3" if ((x >= 23889) & (x <= 23959)) else (
                        "sgmRNA4" if ((x >= 27902) & (x <= 27972)) else (
                            "sgmRNA5" if ((x >= 28285) & (x <= 28355)) else (
                                "sgmRNA6" if ((x >= 28925) & (x <= 28995)) else (
                                    "sgmRNA7" if ((x >= 29622) & (x <= 29692)) else "DVG"
                                )
                            )
                        )
                    )
                ))
            sgmRNAs = forward_junctions[((forward_junctions['start_type'] == "TRSL") & (forward_junctions['stop_type'].str.contains("sgmRNA")))]
            sgmRNA2 = sgmRNAs.loc[sgmRNAs['stop_type'] == "sgmRNA2"].sort_values(by=['depth'], ascending=False)
            sgmRNA3 = sgmRNAs.loc[sgmRNAs['stop_type'] == "sgmRNA3"].sort_values(by=['depth'], ascending=False)
            sgmRNA4 = sgmRNAs.loc[sgmRNAs['stop_type'] == "sgmRNA4"].sort_values(by=['depth'], ascending=False)
            sgmRNA5 = sgmRNAs.loc[sgmRNAs['stop_type'] == "sgmRNA5"].sort_values(by=['depth'], ascending=False)
            sgmRNA6 = sgmRNAs.loc[sgmRNAs['stop_type'] == "sgmRNA6"].sort_values(by=['depth'], ascending=False)
            sgmRNA7 = sgmRNAs.loc[sgmRNAs['stop_type'] == "sgmRNA7"].sort_values(by=['depth'], ascending=False)
            sgmRNA_canonical = pd.DataFrame()
            sgmRNA_canonical = sgmRNA_canonical.append(sgmRNA2.head(1), ignore_index=True)
            sgmRNA_canonical = sgmRNA_canonical.append(sgmRNA3.head(1), ignore_index=True)
            sgmRNA_canonical = sgmRNA_canonical.append(sgmRNA4.head(1), ignore_index=True)
            sgmRNA_canonical = sgmRNA_canonical.append(sgmRNA5.head(1), ignore_index=True)
            sgmRNA_canonical = sgmRNA_canonical.append(sgmRNA6.head(1), ignore_index=True)
            sgmRNA_canonical = sgmRNA_canonical.append(sgmRNA7.head(1), ignore_index=True)
            if (version >= 0.21):
                sgmRNA_alternative = sgmRNAs.merge(sgmRNA_canonical,
                                                   on=['genome', 'start', 'stop', 'type', 'depth', 'strand', 'rgb1', 'rgb2',
                                                       'start_seq', 'stop_seq', 'start_type', 'stop_type', 'frequency',
                                                       'logfreq'], how="left", indicator=True)
            if (version < 0.21):
                sgmRNA_alternative = sgmRNAs.merge(sgmRNA_canonical,
                                                   on=['genome', 'start', 'stop', 'type', 'depth', 'strand', 'start_type', 'stop_type', 'frequency',
                                                       'logfreq'], how="left", indicator=True)
            sgmRNA_alternative = sgmRNA_alternative[sgmRNA_alternative["_merge"] == "left_only"]
            del sgmRNA_alternative["_merge"]
            sgmRNA_alt_summary = pd.DataFrame(columns=['type', 'depth'])
            sgmRNA_alt_summary['type'] = ['sgmRNA2', 'sgmRNA3', 'sgmRNA4', 'sgmRNA5', 'sgmRNA6', 'sgmRNA7']
            sgmRNA2_alt = sgmRNA_alternative.loc[sgmRNA_alternative['stop_type'] == "sgmRNA2"]
            sgmRNA3_alt = sgmRNA_alternative.loc[sgmRNA_alternative['stop_type'] == "sgmRNA3"]
            sgmRNA4_alt = sgmRNA_alternative.loc[sgmRNA_alternative['stop_type'] == "sgmRNA4"]
            sgmRNA5_alt = sgmRNA_alternative.loc[sgmRNA_alternative['stop_type'] == "sgmRNA5"]
            sgmRNA6_alt = sgmRNA_alternative.loc[sgmRNA_alternative['stop_type'] == "sgmRNA6"]
            sgmRNA7_alt = sgmRNA_alternative.loc[sgmRNA_alternative['stop_type'] == "sgmRNA7"]
            sgmRNA2_alt_depth = sum(sgmRNA2_alt['depth'])
            sgmRNA3_alt_depth = sum(sgmRNA3_alt['depth'])
            sgmRNA4_alt_depth = sum(sgmRNA4_alt['depth'])
            sgmRNA5_alt_depth = sum(sgmRNA5_alt['depth'])
            sgmRNA6_alt_depth = sum(sgmRNA6_alt['depth'])
            sgmRNA7_alt_depth = sum(sgmRNA7_alt['depth'])
            sgmRNA_alt_summary.loc[sgmRNA_alt_summary['type'] == "sgmRNA2", ["depth"]] = sgmRNA2_alt_depth
            sgmRNA_alt_summary.loc[sgmRNA_alt_summary['type'] == "sgmRNA3", ["depth"]] = sgmRNA3_alt_depth
            sgmRNA_alt_summary.loc[sgmRNA_alt_summary['type'] == "sgmRNA4", ["depth"]] = sgmRNA4_alt_depth
            sgmRNA_alt_summary.loc[sgmRNA_alt_summary['type'] == "sgmRNA5", ["depth"]] = sgmRNA5_alt_depth
            sgmRNA_alt_summary.loc[sgmRNA_alt_summary['type'] == "sgmRNA6", ["depth"]] = sgmRNA6_alt_depth
            sgmRNA_alt_summary.loc[sgmRNA_alt_summary['type'] == "sgmRNA7", ["depth"]] = sgmRNA7_alt_depth
            DVGs = forward_junctions.loc[((forward_junctions['start_type'] == "TRSL") & (forward_junctions['stop_type'] == "DVG")) | ((forward_junctions['start_type'] == "DVG"))]
            DVGs.to_csv(save_dir_sgmRNAs + sample_name + "_DVGs.txt", sep="\t", index=False)
            sgmRNA_alt_summary.to_csv(save_dir_sgmRNAs + sample_name + "_minor_sgmRNA_summary.txt", sep="\t",
                                      index=False)
            sgmRNA_alternative.to_csv(save_dir_sgmRNAs + sample_name + "_minor_sgmRNAs.txt", sep="\t",
                                      index=False)
            sgmRNA_canonical.to_csv(save_dir_sgmRNAs + sample_name + "_major_sgmRNAs.txt", sep="\t",
                                    index=False)
            sgmRNA_depth = sum(sgmRNAs['depth'])
            sgmRNA_canonical_depth = sum(sgmRNA_canonical['depth'])
            sgmRNA_alternative_depth = sum(sgmRNA_alternative['depth'])
            DVGs_depth = sum(DVGs['depth'])
            sgmRNA_report.loc[sgmRNA_report['sample'] == sample_name, ["total_sgmRNA_depth"]] = sgmRNA_depth
            sgmRNA_report.loc[sgmRNA_report['sample'] == sample_name, ["major_sgmRNA_depth"]] = sgmRNA_canonical_depth
            sgmRNA_report.loc[sgmRNA_report['sample'] == sample_name, ["minor_sgmRNA_depth"]] = sgmRNA_alternative_depth
            sgmRNA_report.loc[sgmRNA_report['sample'] == sample_name, ["DVG_depth"]] = DVGs_depth
        if (virus == "MERS"):
            forward_junctions['start_type'] = forward_junctions['start'].apply(
                lambda x: "TRSL" if ((x >= 32) & (x <= 97)) else "DVG")
            forward_junctions['stop_type'] = forward_junctions['stop'].apply(
                lambda x: "sgmRNA2" if ((x >= 21374) & (x <= 21439)) else (
                    "sgmRNA3" if ((x >=25490) & (x <= 25555)) else (
                        "sgmRNA4" if ((x >= 25812) & (x <= 25877)) else (
                            "sgmRNA5" if ((x >= 26802) & (x <= 26867)) else (
                                "sgmRNA6" if ((x >= 27552) & (x <= 27617)) else (
                                    "sgmRNA7" if ((x >= 27807) & (x <= 27872)) else (
                                        "sgmRNA8" if ((x >= 28514) & (x <= 28579)) else "DVG")))))))
            sgmRNAs = forward_junctions[
                ((forward_junctions['start_type'] == "TRSL") & (forward_junctions['stop_type'].str.contains("sgmRNA")))]
            sgmRNA2 = sgmRNAs.loc[sgmRNAs['stop_type'] == "sgmRNA2"].sort_values(by=['depth'], ascending=False)
            sgmRNA3 = sgmRNAs.loc[sgmRNAs['stop_type'] == "sgmRNA3"].sort_values(by=['depth'], ascending=False)
            sgmRNA4 = sgmRNAs.loc[sgmRNAs['stop_type'] == "sgmRNA4"].sort_values(by=['depth'], ascending=False)
            sgmRNA5 = sgmRNAs.loc[sgmRNAs['stop_type'] == "sgmRNA5"].sort_values(by=['depth'], ascending=False)
            sgmRNA6 = sgmRNAs.loc[sgmRNAs['stop_type'] == "sgmRNA6"].sort_values(by=['depth'], ascending=False)
            sgmRNA7 = sgmRNAs.loc[sgmRNAs['stop_type'] == "sgmRNA7"].sort_values(by=['depth'], ascending=False)
            sgmRNA8 = sgmRNAs.loc[sgmRNAs['stop_type'] == "sgmRNA8"].sort_values(by=['depth'], ascending=False)
            sgmRNA_canonical = pd.DataFrame()
            sgmRNA_canonical = sgmRNA_canonical.append(sgmRNA2.head(1), ignore_index=True)
            sgmRNA_canonical = sgmRNA_canonical.append(sgmRNA3.head(1), ignore_index=True)
            sgmRNA_canonical = sgmRNA_canonical.append(sgmRNA4.head(1), ignore_index=True)
            sgmRNA_canonical = sgmRNA_canonical.append(sgmRNA5.head(1), ignore_index=True)
            sgmRNA_canonical = sgmRNA_canonical.append(sgmRNA6.head(1), ignore_index=True)
            sgmRNA_canonical = sgmRNA_canonical.append(sgmRNA7.head(1), ignore_index=True)
            sgmRNA_canonical = sgmRNA_canonical.append(sgmRNA8.head(1), ignore_index=True)
            if (version >= 0.21):
                sgmRNA_alternative = sgmRNAs.merge(sgmRNA_canonical, on=['genome',
                                                                         'start',
                                                                         'stop',
                                                                         'type',
                                                                         'depth',
                                                                         'strand',
                                                                         'rgb1',
                                                                         'rgb2',
                                                                         'start_seq',
                                                                         'stop_seq',
                                                                         'start_type',
                                                                         'stop_type',
                                                                         'frequency',
                                                                         'logfreq'], how="left", indicator=True)
            if (version < 0.21):
                sgmRNA_alternative = sgmRNAs.merge(sgmRNA_canonical, on=['genome',
                                                                         'start',
                                                                         'stop',
                                                                         'type',
                                                                         'depth',
                                                                         'strand',
                                                                         'start_type',
                                                                         'stop_type',
                                                                         'frequency',
                                                                         'logfreq'], how="left", indicator=True)
            sgmRNA_alternative = sgmRNA_alternative[sgmRNA_alternative["_merge"] == "left_only"]
            del sgmRNA_alternative["_merge"]
            sgmRNA_alt_summary = pd.DataFrame(columns=['type', 'depth'])
            sgmRNA_alt_summary['type'] = ['sgmRNA2', 'sgmRNA3', 'sgmRNA4', 'sgmRNA5', 'sgmRNA6', 'sgmRNA7', 'sgmRNA8']
            sgmRNA2_alt = sgmRNA_alternative.loc[sgmRNA_alternative['stop_type'] == "sgmRNA2"]
            sgmRNA3_alt = sgmRNA_alternative.loc[sgmRNA_alternative['stop_type'] == "sgmRNA3"]
            sgmRNA4_alt = sgmRNA_alternative.loc[sgmRNA_alternative['stop_type'] == "sgmRNA4"]
            sgmRNA5_alt = sgmRNA_alternative.loc[sgmRNA_alternative['stop_type'] == "sgmRNA5"]
            sgmRNA6_alt = sgmRNA_alternative.loc[sgmRNA_alternative['stop_type'] == "sgmRNA6"]
            sgmRNA7_alt = sgmRNA_alternative.loc[sgmRNA_alternative['stop_type'] == "sgmRNA7"]
            sgmRNA8_alt = sgmRNA_alternative.loc[sgmRNA_alternative['stop_type'] == "sgmRNA8"]
            sgmRNA2_alt_depth = sum(sgmRNA2_alt['depth'])
            sgmRNA3_alt_depth = sum(sgmRNA3_alt['depth'])
            sgmRNA4_alt_depth = sum(sgmRNA4_alt['depth'])
            sgmRNA5_alt_depth = sum(sgmRNA5_alt['depth'])
            sgmRNA6_alt_depth = sum(sgmRNA6_alt['depth'])
            sgmRNA7_alt_depth = sum(sgmRNA7_alt['depth'])
            sgmRNA8_alt_depth = sum(sgmRNA8_alt['depth'])
            sgmRNA_alt_summary.loc[sgmRNA_alt_summary['type'] == "sgmRNA2", ["depth"]] = sgmRNA2_alt_depth
            sgmRNA_alt_summary.loc[sgmRNA_alt_summary['type'] == "sgmRNA3", ["depth"]] = sgmRNA3_alt_depth
            sgmRNA_alt_summary.loc[sgmRNA_alt_summary['type'] == "sgmRNA4", ["depth"]] = sgmRNA4_alt_depth
            sgmRNA_alt_summary.loc[sgmRNA_alt_summary['type'] == "sgmRNA5", ["depth"]] = sgmRNA5_alt_depth
            sgmRNA_alt_summary.loc[sgmRNA_alt_summary['type'] == "sgmRNA6", ["depth"]] = sgmRNA6_alt_depth
            sgmRNA_alt_summary.loc[sgmRNA_alt_summary['type'] == "sgmRNA7", ["depth"]] = sgmRNA7_alt_depth
            sgmRNA_alt_summary.loc[sgmRNA_alt_summary['type'] == "sgmRNA8", ["depth"]] = sgmRNA8_alt_depth
            DVGs = forward_junctions.loc[
                ((forward_junctions['start_type'] == "TRSL") & (forward_junctions['stop_type'] == "DVG")) | (
                (forward_junctions['start_type'] == "DVG"))]
            DVGs.to_csv(save_dir_sgmRNAs + sample_name + "_DVGs.txt", sep="\t", index=False)
            sgmRNA_alt_summary.to_csv(save_dir_sgmRNAs + sample_name + "_minor_sgmRNA_summary.txt", sep="\t", index=False)
            sgmRNA_alternative.to_csv(save_dir_sgmRNAs + sample_name + "_minor_sgmRNAs.txt", sep="\t",
                                      index=False)
            sgmRNA_canonical.to_csv(save_dir_sgmRNAs + sample_name + "_major_sgmRNAs.txt", sep="\t",
                                      index=False)
            sgmRNA_depth = sum(sgmRNAs['depth'])
            sgmRNA_canonical_depth = sum(sgmRNA_canonical['depth'])
            sgmRNA_alternative_depth = sum(sgmRNA_alternative['depth'])
            DVGs_depth = sum(DVGs['depth'])
            sgmRNA_report.loc[sgmRNA_report['sample'] == sample_name, ["total_sgmRNA_depth"]] = sgmRNA_depth
            sgmRNA_report.loc[sgmRNA_report['sample'] == sample_name, ["major_sgmRNA_depth"]] = sgmRNA_canonical_depth
            sgmRNA_report.loc[sgmRNA_report['sample'] == sample_name, ["minor_sgmRNA_depth"]] = sgmRNA_alternative_depth
            sgmRNA_report.loc[sgmRNA_report['sample'] == sample_name, ["DVG_depth"]] = DVGs_depth
        if (virus == "SARS2"):
            forward_junctions['start_type'] = forward_junctions['start'].apply(
                lambda x: "TRSL" if ((x >= 40) & (x <= 105)) else "DVG")
            forward_junctions['stop_type'] = forward_junctions['stop'].apply(
                lambda x: "sgmRNA2" if ((x >= 21526) & (x <= 21591)) else (
                        "sgmRNA3" if ((x >= 25355) & (x <= 25420)) else (
                            "sgmRNA4" if ((x >= 26207) & (x <= 26272)) else (
                                "sgmRNA5" if ((x >= 26443) & (x <= 26508)) else(
                                    "sgmRNA6" if ((x >= 27011) & (x <= 27076)) else(
                                        "sgmRNA7" if ((x >= 27358) & (x <= 27423)) else(
                                            "sgmRNA8" if ((x >= 27858) & (x <= 27923)) else(
                                                "sgmRNA9" if ((x >= 28230) & (x <= 28295)) else "DVG"
                                            )
                                        )
                                    )
                                )
                            )
                        )
                    )
                )
            sgmRNAs = forward_junctions[
                # ((forward_junctions['start_type'] == "TRSL") & (forward_junctions['stop_type'].str.contains("sgmRNA")))
                ((forward_junctions['start_type'] == "TRSL") & ((forward_junctions['stop_type'] == "sgmRNA2") |
                                                                (forward_junctions['stop_type'] == "sgmRNA3") |
                                                                (forward_junctions['stop_type'] == "sgmRNA4") |
                                                                (forward_junctions['stop_type'] == "sgmRNA5") |
                                                                (forward_junctions['stop_type'] == "sgmRNA6") |
                                                                (forward_junctions['stop_type'] == "sgmRNA7") |
                                                                (forward_junctions['stop_type'] == "sgmRNA8") |
                                                                (forward_junctions['stop_type'] == "sgmRNA9")
                                                                )
                 )
            ]
            sgmRNAs = sgmRNAs.reset_index(drop=True)
            sgmRNA2 = sgmRNAs.loc[sgmRNAs['stop_type'] == "sgmRNA2"].sort_values(by=['depth'], ascending=False)
            sgmRNA3 = sgmRNAs.loc[sgmRNAs['stop_type'] == "sgmRNA3"].sort_values(by=['depth'], ascending=False)
            sgmRNA4 = sgmRNAs.loc[sgmRNAs['stop_type'] == "sgmRNA4"].sort_values(by=['depth'], ascending=False)
            sgmRNA5 = sgmRNAs.loc[sgmRNAs['stop_type'] == "sgmRNA5"].sort_values(by=['depth'], ascending=False)
            sgmRNA6 = sgmRNAs.loc[sgmRNAs['stop_type'] == "sgmRNA6"].sort_values(by=['depth'], ascending=False)
            sgmRNA7 = sgmRNAs.loc[sgmRNAs['stop_type'] == "sgmRNA7"].sort_values(by=['depth'], ascending=False)
            sgmRNA8 = sgmRNAs.loc[sgmRNAs['stop_type'] == "sgmRNA8"].sort_values(by=['depth'], ascending=False)
            sgmRNA9 = sgmRNAs.loc[sgmRNAs['stop_type'] == "sgmRNA9"].sort_values(by=['depth'], ascending=False)
            sgmRNA_canonical = pd.DataFrame()
            sgmRNA_canonical = sgmRNA_canonical.append(sgmRNA2.head(1), ignore_index=True)
            sgmRNA_canonical = sgmRNA_canonical.append(sgmRNA3.head(1), ignore_index=True)
            sgmRNA_canonical = sgmRNA_canonical.append(sgmRNA4.head(1), ignore_index=True)
            sgmRNA_canonical = sgmRNA_canonical.append(sgmRNA5.head(1), ignore_index=True)
            sgmRNA_canonical = sgmRNA_canonical.append(sgmRNA6.head(1), ignore_index=True)
            sgmRNA_canonical = sgmRNA_canonical.append(sgmRNA7.head(1), ignore_index=True)
            sgmRNA_canonical = sgmRNA_canonical.append(sgmRNA8.head(1), ignore_index=True)
            sgmRNA_canonical = sgmRNA_canonical.append(sgmRNA9.head(1), ignore_index=True)
            if (version >= 0.21):
                sgmRNA_alternative = sgmRNAs.merge(sgmRNA_canonical, on=['genome',
                                                                         'start',
                                                                         'stop',
                                                                         'type',
                                                                         'depth',
                                                                         'strand',
                                                                         'rgb1',
                                                                         'rgb2',
                                                                         'start_seq',
                                                                         'stop_seq',
                                                                         'start_type',
                                                                         'frequency',
                                                                         'stop_type',
                                                                         'logfreq'], how="left", indicator=True)
            if (version < 0.21):
                sgmRNA_alternative = sgmRNAs.merge(sgmRNA_canonical, on=['genome',
                                                                         'start',
                                                                         'stop',
                                                                         'type',
                                                                         'depth',
                                                                         'strand',
                                                                         'start_type',
                                                                         'frequency',
                                                                         'stop_type',
                                                                         'logfreq'], how="left", indicator=True)
            sgmRNA_alternative = sgmRNA_alternative[sgmRNA_alternative["_merge"] == "left_only"]
            sgmRNA_alternative = sgmRNA_alternative.reset_index(drop=True)
            del sgmRNA_alternative["_merge"]
            sgmRNA_alt_summary = pd.DataFrame(columns=['type', 'depth'])
            sgmRNA_alt_summary['type'] = ['sgmRNA2',
                                          'sgmRNA3',
                                          'sgmRNA4',
                                          'sgmRNA5',
                                          'sgmRNA6',
                                          'sgmRNA7',
                                          'sgmRNA8',
                                          'sgmRNA9']
            sgmRNA2_alt = sgmRNA_alternative.loc[sgmRNA_alternative['stop_type'] == "sgmRNA2"]
            sgmRNA3_alt = sgmRNA_alternative.loc[sgmRNA_alternative['stop_type'] == "sgmRNA3"]
            sgmRNA4_alt = sgmRNA_alternative.loc[sgmRNA_alternative['stop_type'] == "sgmRNA4"]
            sgmRNA5_alt = sgmRNA_alternative.loc[sgmRNA_alternative['stop_type'] == "sgmRNA5"]
            sgmRNA6_alt = sgmRNA_alternative.loc[sgmRNA_alternative['stop_type'] == "sgmRNA6"]
            sgmRNA7_alt = sgmRNA_alternative.loc[sgmRNA_alternative['stop_type'] == "sgmRNA7"]
            sgmRNA8_alt = sgmRNA_alternative.loc[sgmRNA_alternative['stop_type'] == "sgmRNA8"]
            sgmRNA9_alt = sgmRNA_alternative.loc[sgmRNA_alternative['stop_type'] == "sgmRNA9"]
            sgmRNA2_alt_depth = sum(sgmRNA2_alt['depth'])
            sgmRNA3_alt_depth = sum(sgmRNA3_alt['depth'])
            sgmRNA4_alt_depth = sum(sgmRNA4_alt['depth'])
            sgmRNA5_alt_depth = sum(sgmRNA5_alt['depth'])
            sgmRNA6_alt_depth = sum(sgmRNA6_alt['depth'])
            sgmRNA7_alt_depth = sum(sgmRNA7_alt['depth'])
            sgmRNA8_alt_depth = sum(sgmRNA8_alt['depth'])
            sgmRNA9_alt_depth = sum(sgmRNA9_alt['depth'])
            sgmRNA_alt_summary.loc[sgmRNA_alt_summary['type'] == "sgmRNA2", ["depth"]] = sgmRNA2_alt_depth
            sgmRNA_alt_summary.loc[sgmRNA_alt_summary['type'] == "sgmRNA3", ["depth"]] = sgmRNA3_alt_depth
            sgmRNA_alt_summary.loc[sgmRNA_alt_summary['type'] == "sgmRNA4", ["depth"]] = sgmRNA4_alt_depth
            sgmRNA_alt_summary.loc[sgmRNA_alt_summary['type'] == "sgmRNA5", ["depth"]] = sgmRNA5_alt_depth
            sgmRNA_alt_summary.loc[sgmRNA_alt_summary['type'] == "sgmRNA6", ["depth"]] = sgmRNA6_alt_depth
            sgmRNA_alt_summary.loc[sgmRNA_alt_summary['type'] == "sgmRNA7", ["depth"]] = sgmRNA7_alt_depth
            sgmRNA_alt_summary.loc[sgmRNA_alt_summary['type'] == "sgmRNA8", ["depth"]] = sgmRNA8_alt_depth
            sgmRNA_alt_summary.loc[sgmRNA_alt_summary['type'] == "sgmRNA9", ["depth"]] = sgmRNA9_alt_depth
            DVGs = forward_junctions.loc[
                ((forward_junctions['start_type'] == "TRSL") & (forward_junctions['stop_type'] == "DVG")) | (
                (forward_junctions['start_type'] == "DVG"))]
            DVGs.to_csv(save_dir_sgmRNAs + sample_name + "_DVGs.txt", sep="\t", index=False)
            sgmRNA_alt_summary.to_csv(save_dir_sgmRNAs + sample_name + "_minor_sgmRNA_summary.txt", sep="\t", index=False)
            sgmRNA_alternative.to_csv(save_dir_sgmRNAs + sample_name + "_minor_sgmRNAs.txt", sep="\t",
                                      index=False)
            sgmRNA_canonical.to_csv(save_dir_sgmRNAs + sample_name + "_major_sgmRNAs.txt", sep="\t",
                                      index=False)
            sgmRNA_depth = sum(sgmRNAs['depth'])
            sgmRNA_canonical_depth = sum(sgmRNA_canonical['depth'])
            sgmRNA_alternative_depth = sum(sgmRNA_alternative['depth'])
            DVGs_depth = sum(DVGs['depth'])
            sgmRNA_report.loc[sgmRNA_report['sample'] == sample_name, ["total_sgmRNA_depth"]] = sgmRNA_depth
            sgmRNA_report.loc[sgmRNA_report['sample'] == sample_name, ["major_sgmRNA_depth"]] = sgmRNA_canonical_depth
            sgmRNA_report.loc[sgmRNA_report['sample'] == sample_name, ["minor_sgmRNA_depth"]] = sgmRNA_alternative_depth
            sgmRNA_report.loc[sgmRNA_report['sample'] == sample_name, ["DVG_depth"]] = DVGs_depth
sgmRNA_report['total_junctions'] = sgmRNA_report["total_sgmRNA_depth"] + sgmRNA_report["DVG_depth"]
sgmRNA_report["percent_DVGs"] = (sgmRNA_report["DVG_depth"] / sgmRNA_report["total_junctions"]) * 100
sgmRNA_report['percent_major_sgmRNA'] = (sgmRNA_report['major_sgmRNA_depth'] / sgmRNA_report['total_junctions']) * 100
sgmRNA_report['percent_minor_sgmRNA'] = (sgmRNA_report['minor_sgmRNA_depth'] / sgmRNA_report['total_junctions']) * 100
sgmRNA_report["DVG_jfreq"] = (sgmRNA_report["DVG_depth"] / sgmRNA_report["total_nts"]) * 1000000
sgmRNA_report['major_sgmRNA_jfreq'] = (sgmRNA_report['major_sgmRNA_depth'] / sgmRNA_report['total_nts']) * 1000000
sgmRNA_report['minor_sgmRNA_jfreq'] = (sgmRNA_report['minor_sgmRNA_depth'] / sgmRNA_report['total_nts']) * 1000000
sgmRNA_report.to_csv(save_dir_sgmRNAs + exp + "_sgmRNA_DVG_report.txt", sep="\t", index=False)