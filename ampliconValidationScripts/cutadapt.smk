
from os import path

config = { 'grid_mult' : 'sbatch --ntasks=1 --mem=16GB --time=24:00:00 --partition=cmbr --cpus-per-task=1',
           'cutadapt' : '~/panfs/miniconda3/bin/cutadapt',
           'threads' : 1
}

rule trim_se:
    """
    Trim single-end reads. Currently this uses cutadapt alone, runs in
    parallel mode, and uses a pre-specified adaptor sequence
    ("CTGTCTCTTATA") filtering reads shorter than 36nt.
    """
    params:
        other='',
        adaptor='CTGTCTCTTATA',
        grid_opts=config['grid_mult'],
        cmd=config['cutadapt'],
    input:
        fq='{samp}.fastq'
    output:
        log='{samp}.cutadapt_log',
        trim='{samp}_trimmed.fq'
    threads:
        config['threads']
    shell:
        '{params.cmd}  {params.other} -a {params.adaptor} ' \
        '-o {output.trim} {input.fq} > {output.log}'
