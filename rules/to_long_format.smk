rule remove_vcf_header:
    input:
        vcf = f"{DATA}/{{sample}}.ready.vcf.gz"
    output:
        vcf_no_header = temp(f"{OUTDIR}/{{sample}}.ready_noheader.vcf")
    threads:
        get_resource('remove_vcf_header', 'threads')
    resources:
        mem=get_resource('remove_vcf_header', 'mem'),
        walltime=get_resource('remove_vcf_header', 'walltime')
    log:
        f"{LOGDIR}/remove_vcf_header/{{sample}}.log"
    shell: "zcat {input.vcf} | sed '/^#/d' > {output.vcf_no_header}"


# TODO ¿Revisar longitud de los sufijos. Cuantos pedazos va a haber para el primer cromosoma?
checkpoint split_chr:
    input:
        vcf_no_header = rules.remove_vcf_header.output.vcf_no_header
    output:
        split_folder = directory(f"{OUTDIR}/split/{{sample}}")
    threads:
        get_resource('split_chr', 'threads')
    resources:
        mem=get_resource('split_chr', 'mem'),
        walltime=get_resource('split_chr', 'walltime')
    params:
        split_size = get_params('split_chr', 'num_lines')
    log:
        f"{LOGDIR}/split_chr/{{sample}}.log"
    shell: 'mkdir {output.split_folder}; split -l {params.num_lines} "{input.vcf_no_header}" "{output.split_folder}/{wildcards.sample}.part-"'


rule paste_header:
    input:
        split_file = f"{OUTDIR}/split/{{sample}}/{{sample}}.part-{{part}}",
    output:
        headed_file = f"{OUTDIR}/headed/{{sample}}/{{sample}}.part-{{part}}"
    threads:
        get_resource('default', 'threads')
    resources:
        mem=get_resource('default', 'mem'),
        walltime=get_resource('default', 'walltime')
    params:
        header = get_params('paste_header', 'header_vcf')
    log:
        f"{LOGDIR}/paste_header/{{sample}}/{{sample}}.part-{{part}}.log"
    shell: 'cat {params.header} {input.split_file} > {output.headed_file}'


rule make_long:
    input:
        headed_file = rules.paste_header.output.headed_file
    output:
        long_format_part = temp(f"{OUTDIR}/long_part/{{sample}}/{{sample}}.part-{{part}}")
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
    shell: 'python3 VCF-Simplify-master/VcfSimplify.py SimplifyVCF -toType table -inVCF {input.headed_file} -outFile {output.long_format_part} -infos {params.infos} -formats {params.formats} -preHeader {params.preHeader} -mode long'


rule filter:
    input:
        long_part = rules.make_long.output.long_format_part
    output:
        filter_part = f"{OUTDIR}/filtered_part/{{sample}}/{{sample}}.part-{{part}}"
    threads:
        get_resource('default', 'threads')
    resources:
        mem=get_resource('default', 'mem'),
        walltime=get_resource('default', 'walltime')
    log:
        f"{LOGDIR}/filtered_part/{{sample}}/{{sample}}.part-{{part}}.log"
    # Set priority = 10 (Default = 0) to this rule and the following ones so that, if 'filter' and 'make_long'
    # are available, an instance of 'filter' is executed first and eliminate a possible long format file 
    priority: 10
    shell: "sed '/\.\/\./d;/0\/0/d' {input.long_part} > {output.filter_part}"


rule remove_tsv_header:
    input:
        filter_part = rules.filter.output.filter_part
    output:
        tsv_no_header = f"{OUTDIR}/filter_no_header/{{sample}}/{{sample}}.part-{{part}}"
    threads:
        get_resource('default', 'threads')
    resources:
        mem=get_resource('default', 'mem'),
        walltime=get_resource('default', 'walltime')
    log:
        f"{LOGDIR}/remove_tsv_header/{{sample}}/{{sample}}.part-{{part}}.log"
    priority: 10
    shell: "tail -n +2 {input.filter_part} > {output.tsv_no_header}"


def aggregate_input(wildcards):
    '''
    aggregate the file names of the random number of files
    generated at the scatter step
    '''
    checkpoint_output = checkpoints.split_chr.get(**wildcards).output[0]
    return expand(f"{OUTDIR}/filter_no_header/{wildcards.sample}/{wildcards.sample}.part-{{i}}",
           i=glob_wildcards(os.path.join(checkpoint_output, f"{wildcards.sample}.part-{{i}}")).i)


rule join_files:
    input:
        aggregate_input
    output:
        combined = f"{OUTDIR}/final/{{sample}}.tsv",
    threads:
        get_resource('default', 'threads')
    resources:
        mem=get_resource('default', 'mem'),
        walltime=get_resource('default', 'walltime')
    params:
        header = get_params('join_files', 'header_tsv')
    log:
        f"{LOGDIR}/join_files/{{sample}}.log"
    priority: 10
    shell: "cat {params.header} {input} > {output.combined}"

