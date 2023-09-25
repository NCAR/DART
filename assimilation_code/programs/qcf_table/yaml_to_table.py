import yaml

#Prompt user for name of input and output files
input_yaml = input('Please enter the name of your input yaml file (filename must end in ".yaml") OR press enter/return to use the default filename "qcf_table.yaml"\n')
output_txt = input('Please enter the name for the output file for the table (filename must end in ".txt") OR press enter/return to use the default filename "qcf_table.txt"\n')

#Using deault names for input/output files if not specified
if (input_yaml == ''):
    input_yaml = 'qcf_table.yaml'

if (output_txt == ''):
    output_txt = 'qcf_table.txt'

#Open and load yaml file
with open(input_yaml) as file:
    dict = yaml.safe_load(file)

    column_headers = list(dict.keys())
    column_data = list(dict.values())

    obs_errror_info_header = dict['QTY_TEMPLATE']['obs_error_info']
    probit_inflation_header = dict['QTY_TEMPLATE']['probit_inflation']
    probit_state_header = dict['QTY_TEMPLATE']['probit_state']
    probit_extended_state_header = dict['QTY_TEMPLATE']['probit_extended_state']
    obs_inc_info_header = dict['QTY_TEMPLATE']['obs_inc_info']

    f = open(output_txt, "w")

#Write the table's headers to the output file
    f.write(column_headers[0] + ": " + str(column_data[0]) + "\n")

    f.write(column_headers[1] + ": ")
    for name in obs_errror_info_header:
        f.write(name)
    f.write(" ")
    for name in probit_inflation_header:
        f.write(name)
    f.write(" ")
    for name in probit_state_header:
        f.write(name)
    f.write(" ")
    for name in probit_extended_state_header:
        f.write(name)
    f.write(" ")
    for name in obs_inc_info_header:
        f.write(name)
    f.write("\n")

#Writing table data to the output file 
    for i in range(2, len(column_headers)):
        f.write(column_headers[i] + " ")

        obs_error_info = dict[column_headers[i]]['obs_error_info'].items()
        for key, value in obs_error_info:
            f.write(str(value) + " ")

        probit_inflation = dict[column_headers[i]]['probit_inflation'].items()
        for key, value in probit_inflation:
            f.write(str(value) + " ")

        probit_state = dict[column_headers[i]]['probit_state'].items()
        for key, value in probit_state:
            f.write(str(value) + " ")

        probit_extended_state = dict[column_headers[i]]['probit_extended_state'].items()
        for key, value in probit_extended_state:
            f.write(str(value) + " ")

        obs_inc_info = dict[column_headers[i]]['obs_inc_info'].items()
        for key, value in obs_inc_info:
            f.write(str(value) + " ")

        f.write("\n")

    f.close

print('QCF table produced in ' + output_txt)
