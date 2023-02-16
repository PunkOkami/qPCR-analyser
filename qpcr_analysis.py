
# Main part of program to analyse qPCr results from tsv file
# Copyright (C) 2023  Maks Chmielewski
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import csv
import argparse
from sys import exit as s_exit

parser = argparse.ArgumentParser(prog='qPCR Analysis Tool', description='This is a simple command line tool for basic analysis of qPCR results')
parser.add_argument('file_name', help='path to the file with qpcr run results in tsv format')
parser.add_argument('control_sample_name', help='name of control sample in data file')
parser.add_argument('stress_sample_name', help='name of stress sample in data file')
parser.add_argument('tested_genes', help='comma seperated list of tested genes Eg: "gene1,gene2,gene3"')
parser.add_argument('housekeeping_gene', help='housekeeping gene to be used for normalization of ct values')
parser.add_argument('prim_3_quality_gene', help='3 prim part of gene used for quality control of samples')
parser.add_argument('prim_5_quality_gene', help='5 prim part of gene used for quality control of samples')
args = parser.parse_args()

file = open(args.file_name)
reader = csv.reader(file, delimiter='\t')

control_sample_name = args.control_sample_name
stress_sample_name = args.stress_sample_name
tested_genes = args.tested_genes.split(',')
housekeeping_gene = args.housekeeping_gene
prim_3_quality_gene = args.prim_3_quality_gene
prim_5_quality_gene = args.prim_5_quality_gene

sample_name_index = 0
detector_name_index = 0
ct_value_index = 0
control_data = {}
stress_data = {}

for row in reader:
	if len(row) < 2:
		continue
	if 'Well' in row:
		sample_name_index = row.index('Sample Name')
		detector_name_index = row.index('Detector Name')
		ct_value_index = row.index('Ct')
	if row[sample_name_index] == control_sample_name:
		gene_name = row[detector_name_index]
		ct_value = row[ct_value_index]
		if ct_value == 'Undetermined':
			ct_value = -1.0
		else:
			ct_value = float(ct_value)
		ct_values = control_data.get(gene_name, [])
		ct_values.append(ct_value)
		control_data[gene_name] = ct_values
	if row[sample_name_index] == stress_sample_name:
		gene_name = row[detector_name_index]
		ct_value = row[ct_value_index]
		if ct_value == 'Undetermined':
			ct_value = -1.0
		else:
			ct_value = float(ct_value)
		ct_values = stress_data.get(gene_name, [])
		ct_values.append(ct_value)
		stress_data[gene_name] = ct_values

for (gene_name, ct_values) in control_data.items():
	ct_values = [ct for ct in ct_values if ct != -1.0]
	ct_values = round(sum(ct_values)/len(ct_values), 7)
	control_data[gene_name] = ct_values

for (gene_name, ct_values) in stress_data.items():
	ct_values = [ct for ct in ct_values if ct != -1.0]
	ct_values = round(sum(ct_values)/len(ct_values), 7)
	stress_data[gene_name] = ct_values


print(f'Gene used as referance is {housekeeping_gene}\n')
delta_ct_values_control = {}
if len(control_data) == 0:
	print('FATAL ERROR:\nControl sample is not available, calculation of control delta Ct halted.')
	print('Check the name of control sample with data file')
	s_exit(1)
elif control_data.get(housekeeping_gene, 'Not there') == 'Not there':
	print('FATAL ERROR:\nThere is no Ct data available for housekeeping gene. Execution halted. Data for other genes available.')
	print('Check if the name of housekeeping gene is correct')
	s_exit(1)
elif control_data[housekeeping_gene] == -1.0:
	print('FATAL ERROR:\nCt value of houseeekping gene is "Undetermined". qPCR machine could not detect any amplification of housekeeping gene in control samples.')
	print('Redo the experiment')
	s_exit(1)
else:
	print('Control data:')
	ct_prim_3 = control_data.get(prim_3_quality_gene, 'Not there')
	ct_prim_5 = control_data.get(prim_5_quality_gene, 'Not there')
	if ct_prim_3 == 'Not there' or ct_prim_5 == "Not there":
		print('Name of one of quality genes is not available, compare the names with data file')
	elif ct_prim_3 == -1.0 or ct_prim_5 == -1.0:
		print('Ct values for at least one of quality genes is "Undetermined". Quality ratio cannot be calculated.')
		print('Redoing the experiment recommended')
	else:
		quality_ratio = control_data[prim_3_quality_gene]/control_data[prim_5_quality_gene]
		print(f'Quality ratio = {quality_ratio}')
	print('Delta Ct values for each of tested genes:\nGene --- Ct value --- Delta Ct value')
	for gene in tested_genes:
		control_ct_gene = control_data.get(gene, 'Not there')
		if control_ct_gene == 'Not there':
			print(f'Gene name {gene} was not found in control sample data, check if input name is correct')
		elif control_ct_gene == -1.0:
			print(f'Ct values for {gene} are "Undetermined". No amplification of this gene was detected. Delta Ct cannot be calculated')
			delta_ct = 'Undetermined'
			delta_ct_values_control[gene] = delta_ct
		else:
			delta_ct = control_ct_gene - control_data[housekeeping_gene]
			print(f'{gene} --- {control_ct_gene} -- {delta_ct}')
			delta_ct_values_control[gene] = delta_ct

print('\n')
if len(delta_ct_values_control) == 0:
	print('FATAL ERROR:\nThere are no Delta Ct values for genes in control sample data. No genes in tested genes were found in data file')
	s_exit(1)

control_usable_data = 0
for delta_ct in delta_ct_values_control.values():
	if delta_ct != "Undetermined":
		control_usable_data += 1
if control_usable_data == 0:
	print('FATAL ERROR:\nNo genes in control sample had Ct values. Calculation of Delta Ct failed, program halted')
	s_exit(1)
else:
	print(f'{control_usable_data} out of {len(tested_genes)} input genes in control sample have Delta Ct calculated.', end='')
	print(f'Those genes will be used for calculation of Delta Delta Ct values')
	print('\n')

delta_ct_values_stress = {}
if len(stress_data) == 0:
	print('FATAL ERROR:\nStress sample is not available, calculation of stress delta Ct halted.')
	print('Check the name of control sample with data file')
	s_exit(1)
elif stress_data.get(housekeeping_gene, 'Not there') == 'Not there':
	print('FATAL ERROR:\nThere is no Ct data available for housekeeping gene. Execution halted. Data for other genes available.')
	print('Check if the name of housekeeping gene is correct')
	s_exit(1)
elif stress_data[housekeeping_gene] == -1.0:
	print('FATAL ERROR:\nCt value of houseeekping gene is "Undetermined". qPCR machine could not detect any amplification of housekeeping gene in stress samples.')
	print('Redo the experiment')
	s_exit(1)
else:
	print('Stress data:')
	ct_prim_3 = stress_data.get(prim_3_quality_gene, 'Not there')
	ct_prim_5 = stress_data.get(prim_5_quality_gene, 'Not there')
	if ct_prim_3 == 'Not there' or ct_prim_5 == "Not there":
		print('Name of one of quality genes is not available, compare the names with data file')
	elif ct_prim_3 == -1.0 or ct_prim_5 == -1.0:
		print('Ct values for at least one of quality genes is "Undetermined". Quality ratio cannot be calculated.')
		print('Redoing the experiment recommended')
	else:
		quality_ratio = stress_data[prim_3_quality_gene]/stress_data[prim_5_quality_gene]
		print(f'Quality ratio = {quality_ratio}')
	print('Delta Ct values for each of tested genes:\nGene --- Ct value --- Delta Ct value')
	for gene in tested_genes:
		stress_ct_gene = stress_data.get(gene, 'Not there')
		if stress_ct_gene == 'Not there':
			print(f'Gene name {gene} was not found in control sample data, check if input name is correct')
		elif stress_ct_gene == -1.0:
			print(f'Ct values for {gene} are "Undetermined". No amplification of this gene was detected. Delta Ct cannot be calculated')
			delta_ct = 'Undetermined'
			delta_ct_values_stress[gene] = delta_ct
		else:
			delta_ct = stress_ct_gene - stress_data[housekeeping_gene]
			print(f'{gene} --- {stress_ct_gene} -- {delta_ct}')
			delta_ct_values_stress[gene] = delta_ct

print('\n')
if len(delta_ct_values_stress) == 0:
	print('FATAL ERROR:\nThere are no Delta Ct values for genes in stress sample data. No genes in tested genes were found in data file')
	s_exit(1)

stress_usable_data = 0
for delta_ct in delta_ct_values_stress.values():
	if delta_ct != "Undetermined":
		stress_usable_data += 1
if stress_usable_data == 0:
	print('FATAL ERROR:\nNo genes in stress sample had Ct values. Calculation of Delta Ct failed, program halted')
	s_exit(1)
else:
	print(f'{control_usable_data} out of {len(tested_genes)} input genes in control sample have Delta Ct calculated.', end='')
	print(f'Those genes will be used for calculation of Delta Delta Ct values')
	print('\n')

print('Delta Delta Ct values for each of tested genes:')
delta_delta_ct_values = {}
for gene in tested_genes:
	delta_ct_control = delta_ct_values_control.get(gene, 'Not there')
	delta_ct_stress = delta_ct_values_stress.get(gene, 'Not there')
	if delta_ct_control in ('Undetermined', 'Not there') or delta_ct_stress in ('Undetermined', 'Not there'):
		print(f'Gene {gene} had incorrect Ct value in at least one of samples, check output above for details')
	else:
		delta_delta_ct = delta_ct_stress - delta_ct_control
		delta_delta_ct_values[gene] = delta_delta_ct
		print(f'{gene} --- {delta_delta_ct}')

print('\n')
if len(delta_delta_ct_values) == 0:
	print('FATAL ERROR\nThere are no calculated Delta Delta Ct values, check the output above for details')
	s_exit(1)

print('Fold change for each of tested genes:')
for (gene, delta_delta_ct) in delta_delta_ct_values.items():
	delta_delta_ct = -delta_delta_ct
	fold_change = 2**delta_delta_ct
	print(f'{gene} --- {fold_change}')

increased_genes = [key for (key, value) in delta_delta_ct_values.items() if value < 0]
decreased_genes = [key for (key, value) in delta_delta_ct_values.items() if value > 0]

print('\n')
if increased_genes:
	print("These genes had bigger expression in stress samples, which suggests they take part in organism's response to this particular stress")
	for gene in increased_genes:
		print(f' - {gene}')

print('\n')
if decreased_genes:
	print("These genes had lower expression in stress samples, which suggests they would hamper organism's efforts in responding to this particular stress")
	for gene in decreased_genes:
		print(f' - {gene}')
