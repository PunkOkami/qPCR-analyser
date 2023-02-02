
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

tested_genes = ['ATG8H', 'HSP101', 'RCAR3']
housekeeping_gene = 'ef1alfa'
prim_3_quality_gene = 'GAPDH3'
prim_5_quality_gene = 'GAPDH5'
print(f'Gene used as referance is {housekeeping_gene}\n')

delta_ct_values_control = {}
if len(control_data) == 0:
	print('FATAL ERROR:\nControl sample is not available, calculation of control delta Ct halted')
elif control_data.get(housekeeping_gene, 'Not there') == 'Not there':
	print('FATAL ERROR:\nThere is no Ct data available for housekeeping gene, execution halted, but there is some data for other genes.')
	print('Check if name of housekeeping gene is correct')
elif control_data[housekeeping_gene] == -1.0:
	print('FATAL ERROR:\nCt value is Undetermined. qPCR machine could not see any amplification of housekeeping gene in control samples.')
	print('Redo the experiment')
else:
	print('Control data:')
	ct_prim_3 = control_data.get(prim_3_quality_gene, 'Not there')
	ct_prim_5 = control_data.get(prim_5_quality_gene, 'Not there')
	if ct_prim_3 == 'Not there' or ct_prim_5 == "Not there":
		print('Name of one of quality genes is not available, check the names with data file')
	elif ct_prim_3 == -1.0 or ct_prim_5 == -1.0:
		print('Ct values for one of quality genes are Undetermined. Quality ratio cannot be calculated, you might want to redo the experiment')
	else:
		quality_ratio = control_data[prim_3_quality_gene]/control_data[prim_5_quality_gene]
		print(f'Quality ratio = {quality_ratio}')
	print('Delta Ct values for each of tested genes:\nGene --- Ct value --- Delta Ct value')
	delta_ct = 0
	for gene in tested_genes:
		control_ct_gene = control_data.get(gene, 'Not there')
		if control_ct_gene == 'Not there':
			print(f'Gene name {gene} was not found in data for control samples, check if name is correct')
		elif control_ct_gene == -1.0:
			print(f"Ct values for {gene} are Undetermined, there was no amplification of this gene, delta_ct cannot be calculated")
			delta_ct = 'Undetermined'
		else:
			delta_ct = control_ct_gene - control_data[housekeeping_gene]
			print(f'{gene} --- {control_ct_gene} -- {delta_ct}')
		delta_ct_values_control[gene] = delta_ct

print('\n')
delta_ct_values_stress = {}
print('Stress data:')
if stress_data[prim_3_quality_gene] == -1.0 or stress_data[prim_5_quality_gene] == -1.0:
	print('Ct values for at least one of quality genes are not avaible. Quality ratio cannot be calculated, consider redoing the epreriment')
else:
	quality_ratio = stress_data[prim_3_quality_gene]/stress_data[prim_5_quality_gene]
	print(f'Quality ratio = {quality_ratio}')
print('Delta Ct values for each of tested genes:')
for gene in tested_genes:
	delta_ct = stress_data[gene] - stress_data[housekeeping_gene]
	delta_ct_values_stress[gene] = delta_ct
	print(f'{gene} --- {delta_ct}')

print('\n')
print('Delta delta Ct values for each tested genes:')
delta_delta_ct_values = {}
for gene in tested_genes:
	delta_delta_ct = delta_ct_values_stress[gene] - delta_ct_values_control[gene]
	delta_delta_ct_values[gene] = delta_delta_ct
	print(f'{gene} --- {delta_delta_ct}')


print('\n')
print('Fold change for each of tested genes:')
for gene in tested_genes:
	delta_delta_ct = -delta_delta_ct_values[gene]
	fold_change = 2**delta_delta_ct
	print(f'{gene} --- {fold_change}')

increased_genes = [key for (key, value) in delta_delta_ct_values.items() if value < 0]
decreased_genes = [key for (key, value) in delta_delta_ct_values.items() if value > 0]

print('\n')
if increased_genes:
	print("These genes had bigger expression in stress samples, which means they take part in organism's response to this particular stress")
	for gene in increased_genes:
		print(f' - {gene}')

print('\n')
if decreased_genes:
	print("These genes had lower expression in stress samples, which means they would hamper organism's efforts in responding to this particular stress")
	for gene in decreased_genes:
		print(f' - {gene}')
