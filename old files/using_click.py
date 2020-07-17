
#!/usr/bin/env python

import os
import subprocess
import click

#python using_fastani.py --inputs /mnt/c/Users/Eilish\ McMaster/PycharmProjects/kea/test_files --output output_dir fastani
#python using_fastani.py fastani --inputs ~/kea_wd/test_files --output output_dir

@click.group()
@click.option('--output', type=click.Path(), help='Directory for output polished MAG')
@click.option('--inputs', type=click.Path(exists=True), help='Absolute path of directory containing MAGs')

#@click.option('--threads', default=16, help='Number of threads to use')

#make output directory
def setup(output, inputs):
    #make output directory as specified by --output
    os.makedirs(output)
    os.chdir(output)
    #make a working directory
    os.makedirs('kea_wd')
    os.chdir('kea_wd')

    #make a list of all of the input files with their absolute path to use later
    files_in_input_dir = os.listdir(inputs)
    input_mag_abs_path = open('input_mag_abs_path.tsv', 'w')
    for input_file in files_in_input_dir:
        input_mag_abs_path.write("%s/%s\n" % (inputs, input_file))
    input_mag_abs_path.close()

@setup.command()
def fastani():
    click.echo("Is fastANI running?")
    subprocess.Popen("fastANI --ql input_mag_abs_path.tsv --rl input_mag_abs_path.tsv -t 20 -o ani_output.tsv" , shell=True).wait()

if __name__ == "__main__":
    setup()
    fastani()