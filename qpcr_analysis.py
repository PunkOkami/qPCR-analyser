import csv

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
		elif  ct_mean == 'Undertermined':
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
print('Control data:')
if control_data[prim_3_quality_gene] == -1.0 or control_data[prim_5_quality_gene] == (-1.0):
	print('Ct values for one of quality genes are not avaible. Quality ratio cannot be calculated, consider redoing the epreriment')
else:
	quality_ratio = control_data[prim_3_quality_gene]/control_data[prim_5_quality_gene]
	print(f'Quality ratio = {quality_ratio}')
print('Delta Ct values for each of tested genes:')
for gene in tested_genes:
	delta_ct = control_data[gene] - control_data[housekeeping_gene]
	delta_ct_values_control[gene] = delta_ct
	print(f'{gene} --- {delta_ct}')

print('\n')
delta_ct_values_stress = {}
print('Stress data:')
quality_ratio = float(stress_data['GAPDH3'])/float(stress_data['GAPDH5'])
print(f'Quality ratio = {quality_ratio}')
print('Delta Ct values for each of tested genes:')
for gene in tested_genes:
	delta_ct = float(stress_data[gene]) - float(stress_data[housekeeping_gene])
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


# ToDo:
# - make chosing fields in tsv file based on names, not idexes
# - add command line options
# - add testing for empty values in tsv file
# - pack delta Ct calculation into a function, less code

