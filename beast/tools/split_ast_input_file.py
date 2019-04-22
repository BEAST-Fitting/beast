import numpy as np
import subprocess

def split_asts(project, giant_ast_file, n_per_file):
    """
    Divide the list of ASTs into sub-files, to make it easier to process them
    with Ben's pipeline.  The sub-files will be placed in the project folder,
    and also put into a tar file.

    Parameters
    ----------
    project : string
        name of the project from datamodel.project, which defines the folder name

    giant_ast_file : string
        Name of the AST file created by AST section of run_beast.py

    n_per_file : integer
        number of lines per output file

    """

    with open(giant_ast_file,'r') as f:

        ast_data = f.readlines()

        # remove the first line
        del ast_data[0]
 
        # length of file
        n_lines = len(ast_data)
        print('AST file contains '+str(n_lines)+' lines')
    
        # number of new files
        n_file = n_lines // n_per_file
            
        # check if there are extras
        if n_lines % n_per_file != 0:
            n_file += 1

        print('Splitting AST file into '+str(n_file)+' sub-files')

        # loop across files
        for i in range(n_file):

            #print('writing file ', i+1, ' of ', n_file)

            with open('./'+project+'/fake_'+str(i+1)+'.lst','w') as f_out:

                # write n_per_file lines from large file
                for j in range(i*n_per_file, (i+1)*n_per_file):

                    # make sure we don't go over for the last file
                    if j > n_lines-1:
                        break
                
                    f_out.write(ast_data[j])
            


    # combine AST files into a tar file

    file_list = ['./'+project+'/fake_'+str(i+1)+'.lst' for i in range(n_file)]

    cmd = 'tar -cf ' + './'+project+'/'+project + '_ASTs.tar ' + ' '.join(file_list)

    subprocess.run(cmd, shell=True)

    # remove the individual files to avoid clutter

    cmd = 'rm '+' '.join(file_list)
    subprocess.run(cmd, shell=True)
