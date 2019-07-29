
from os import path

config = {'grid_single' : 'sbatch --ntasks=1 --mem=16GB --time=00:30:00 --partition=cmbr --cpus-per-task=1'}

rule filter:
    """
    Run filtering on the meth data
    """
    params:
        grid_opts=config['grid_single']
    input:
        meth='{samp}.meth'
    output:
        meth='{samp}.meth.highcoverage'
    shell:
        'awk -F "\t" \'$6 > 1000\' {input.meth} > {output.meth}'
