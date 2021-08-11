rule remove_header:
    input:
        vcf = f"{DATA}/{{sample}}.ready.vcf.gz"
    output:
        vcf_no_header = temp(f"{OUTDIR}/{{sample}}.ready_noheader.vcf")
    threads:
        get_resource('default', 'threads')
    resources:
        mem=get_resource('default', 'mem'),
        walltime=get_resource('default', 'walltime')
    log:
        f"{LOGDIR}/remove_header/{{sample}}.log"
    shell: "zcat {input.vcf} | sed '/^#/d' > {output.vcf_no_header}"


# TODO ¿Revisar longitud de los sufijos. Cuantos pedazos va a haber para el primer cromosoma?
checkpoint split_chr:
    input:
        vcf_no_header = rules.remove_header.output.vcf_no_header
    output:
        split_folder = directory(f"{OUTDIR}/split/{{sample}}")
    threads:
        get_resource('default', 'threads')
    resources:
        mem=get_resource('default', 'mem'),
        walltime=get_resource('default', 'walltime')
    log:
        f"{LOGDIR}/split_chr/{{sample}}.log"
    shell: 'mkdir {output.split_folder}; split -d -l 50000 "{input.vcf_no_header}" "{output.split_folder}/{wildcards.sample}.part-"'


rule paste_header:
    input:
        #split_folder = rules.split_chr.output.split_folder,
        split_file = f"{OUTDIR}/split/{{sample}}/{{sample}}.part-{{part}}",
        header = get_params('paste_header', 'header')
    output:
        headed_file = f"{OUTDIR}/headed/{{sample}}/{{sample}}.part-{{part}}"
    threads:
        get_resource('default', 'threads')
    resources:
        mem=get_resource('default', 'mem'),
        walltime=get_resource('default', 'walltime')
    log:
        f"{LOGDIR}/paste_header/{{sample}}/{{sample}}.part-{{part}}.log"
    shell: 'cat {input.header} {input.split_file} > {output.headed_file}'


rule make_long:
    input:
        headed_file = rules.paste_header.output.headed_file
    output:
        long_format_part = f"{OUTDIR}/long_part/{{sample}}/{{sample}}.part-{{part}}"
    threads:
        get_resource('make_long', 'threads')
    resources:
        mem=get_resource('make_long', 'mem'),
        walltime=get_resource('make_long', 'walltime')
    params:
        infos = get_params('make_long', 'infos'),
        formats = get_params('make_long', 'formats'),
        preHeader = get_params('make_long', 'preHeader')
    log:
        f"{LOGDIR}/make_long/{{sample}}/{{sample}}.part-{{part}}.log"
    shell: 'python3 VCF-Simplify/VcfSimplify.py SimplifyVCF -toType table -inVCF {input.headed_file} -outFile {output.long_format_part} -infos {params.infos} -formats {params.formats} -preHeader {params.preHeader} -mode long'


def aggregate_input(wildcards):
    '''
    aggregate the file names of the random number of files
    generated at the scatter step
    '''
    checkpoint_output = checkpoints.split_chr.get(**wildcards).output[0]
    aux = expand(f"{OUTDIR}/long_part/{wildcards.sample}/{wildcards.sample}.part-{{i}}",
           i=glob_wildcards(os.path.join(checkpoint_output, f"{wildcards.sample}.part-{{i}}")).i)
    print(aux)
    return aux


rule join_files:
    input:
        aggregate_input
    output:
        combined = f"{OUTDIR}/final/{{sample}}.tsv",
    shell:
        'cat {input} > {output.combined}'

