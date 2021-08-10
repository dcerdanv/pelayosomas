
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
rule split_chr:
    input:
        vcf_no_header = rules.remove_header.output.vcf_no_header
    output:
        split_folder=directory(f"{OUTDIR}/split/{{sample}}")
    threads:
        get_resource('default', 'threads')
    resources:
        mem=get_resource('default', 'mem'),
        walltime=get_resource('default', 'walltime')
    log:
        f"{LOGDIR}/split_chr/{{sample}}.log"
    shell: 'mkdir {output.split_folder}; split -l 50000 "{input.vcf_no_header}" "{output.split_folder}/{wildcards.sample}.part-"'

