
from os import path

config = {'grid_single' : 'sbatch --ntasks=1 --mem=16GB --time=01:00:00 --partition=cmbr --cpus-per-task=1'}

rule methcounts:
    """
    Run methcounts on the amplicon data
    """
    params:
        grid_opts=config['grid_single']
    input:
        mr='{samp}.mr.sorted'
    output:
        meth='{samp}.meth'
    shell:
        '/home/rcf-47/andrewds/as/code/for_peter/methpipe/bin/methcounts -v -c ~/staging/chroms -n -o {output.meth} {input.mr}'
