from os import path

config = {'grid_single' : 'sbatch --ntasks=1 --mem=16GB --time=02:00:00 --partition=cmbr --cpus-per-task=1'}

rule filter:
    """
    Run filtering on the meth data
    """
    params:
        grid_opts=config['grid_single']
    input:
        meth='{samp}.meth'
    output:
        meth='{samp}.meth.roi'
    shell:
        '/home/rcf-47/andrewds/as/code/for_peter/methpipe/bin/selectsites -v 13_sites.bed -o {output.meth} {input.meth}'


