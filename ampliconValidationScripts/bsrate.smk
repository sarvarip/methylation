from os import path

config = {'grid_single' : 'sbatch --ntasks=1 --mem=16GB --time=02:00:00 --partition=cmbr --cpus-per-task=1'}

rule filter:
    """
    Run filtering on the meth data
    """
    params:
        grid_opts=config['grid_single']
    input:
        mr='{samp}.mr.sorted'
    output:
        bsrate='{samp}.bsrate'
    shell:
        '/home/rcf-47/andrewds/as/code/for_peter/methpipe/bin/bsrate -c ~/staging/chroms -o {output.bsrate} {input.mr}'


