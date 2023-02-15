
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
from sys import exit as s_exit

file = open('data/qPCR_results_custom_genes.txt')
reader = csv.reader(file, delimiter = '\t')

sample_name_index = 1
detector_name_index = 2
ct_mean_index = 12
control_data = {}
stress_data = {}
for row in reader:
	if len(row) < 2:
		continue
	if row[1] == 'kontrola Tomek i Zuzia':
		gene_name = row[detector_name_index]
		ct_mean = row[ct_mean_index]
		if ct_mean == '':
			ct_mean = -1.0
		elif  ct_mean == 'Undertermined':
			ct_mean = -1.0
		else:
			ct_mean = float(ct_mean)
		control_data[gene_name] = ct_mean
	if row[1] == 'stres Tomek i Zuzio':
		gene_name = row[detector_name_index]
		ct_mean = row[ct_mean_index]
		if ct_mean == '':
			ct_mean = -1.0
		else:
			ct_mean = float(ct_mean)
		stress_data[gene_name] = ct_mean

tested_genes = ['ATG8h', 'HSP101', 'RCAR3']
housekeeping_gene = 'ef1alfa'
prim_3_quality_gene = 'GAPDH3'
prim_5_quality_gene = 'GAPDH5'
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
