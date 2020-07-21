import os
def input_click(input, output):
    # Make output directory as specified by --output
    if os.path.exists(output + '/kea_wd'):
        os.chdir(output + '/kea_wd')
        print('Output path already exists')
    else:
        os.makedirs(output)
        os.chdir(output)
        print('Output directory created')
        # make a working directory
        os.makedirs('kea_wd')
        os.chdir('kea_wd')
        print('Working directory created')

    # Make a list of all of the input files with their absolute path to use later
    if os.path.isfile('input_mag_abs_path.tsv'):
        print('Input path files already exists')
    else:
        files_in_input_dir = os.listdir(input)
        input_mag_abs_path = open('start.tsv', 'w')
        for input_file in files_in_input_dir:
            input_mag_abs_path.write("%s/%s\n" % (input, input_file))
        input_mag_abs_path.close()
        # remove spaces
        fin = open('start.tsv', 'rt')
        fout = open('input_mag_abs_path.tsv', 'wt')
        for line in fin:
            fout.write(line.replace(' ', '\ '))
        fin.close()
        fout.close()
        print('File containing list of input files with absolute path created (input_mag_abs_path.tsv)')