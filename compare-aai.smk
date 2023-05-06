import pandas as pd

out_dir = "output.pocp"
logs_dir = os.path.join(out_dir, 'logs')
basename = 'brady'
ksizes = [6,7,8,9,10,11]
scaled = [5,10,20,40,100]
moltype='protein'
POCP_TABLE="brady_pocp_table.tab"

rule all:
    input: 
        expand(os.path.join(out_dir, f"{basename}.{{moltype}}.sc{{scaled}}.zip"), moltype=moltype, scaled=scaled),
        expand(os.path.join(out_dir, f"{basename}.{{moltype}}.k{{k}}-sc{{scaled}}.compare.csv"), moltype=moltype, scaled=scaled, k=ksizes),
        expand(os.path.join(out_dir, "plots", f"{basename}.{{moltype}}.k{{k}}-sc{{scaled}}.compare.matrix.png"), moltype=moltype, scaled=scaled, k=ksizes),
        expand(os.path.join(out_dir, "pocp-compare", f"{basename}.{{moltype}}.k{{k}}-sc{{scaled}}.pocp-compare.csv"), moltype=moltype, scaled=scaled, k=ksizes),

def make_param_str(ksizes, scaled):
    ks = [ f'k={k}' for k in ksizes ]
    ks = ",".join(ks)
    if isinstance(scaled, list):
        scaled = min(scaled) #take minimum value of scaled list
    return f"{ks},scaled={scaled},abund"


rule sketch_fromfile:
    input: 
        fromfile=os.path.join(out_dir, "{basename}.fromfile.csv"),
    output: os.path.join(out_dir, "{basename}.{moltype}.zip")
    params:
        sketch_params=make_param_str(ksizes=ksizes, scaled=scaled),
    threads: 1
    resources:
        mem_mb=6000,
        time=4000,
        partition="high2",
    conda: "conf/env/sourmash.yml"
    log:  os.path.join(logs_dir, "sketch", "{basename}.{moltype}.log")
    benchmark:  os.path.join(logs_dir, "sketch", "{basename}.{moltype}.benchmark")
    shell:
        """
        sourmash sketch fromfile {input.fromfile} -p {params.sketch_params} -o {output} 2> {log}
        """

rule downsample:
    input: 
        zip=os.path.join(out_dir, "{basename}.{moltype}.zip")
    output: 
        zip=os.path.join(out_dir, "{basename}.{moltype}.sc{scaled}.zip")
    params:
        alpha_cmd=lambda w: f"--{w.moltype}"
    threads: 1
    resources:
        mem_mb=6000,
        time=4000,
        partition="high2",
    conda: "conf/env/sourmash.yml"
    log:  os.path.join(logs_dir, "downsample", "{basename}.{moltype}.sc{scaled}.log")
    benchmark:  os.path.join(logs_dir, "downsample", "{basename}.{moltype}.sc{scaled}.benchmark")
    shell:
        """
        sourmash sig downsample {input.zip} --scaled {wildcards.scaled} {params.alpha_cmd} -o {output} 2> {log}
        """


rule compare_aai:
    input: 
        zip=os.path.join(out_dir, "{basename}.{moltype}.sc{scaled}.zip")
    output: 
        csv=os.path.join(out_dir, "{basename}.{moltype}.k{k}-sc{scaled}.compare.csv"),
        np=os.path.join(out_dir, "{basename}.{moltype}.k{k}-sc{scaled}.compare"),
        labels=os.path.join(out_dir, "{basename}.{moltype}.k{k}-sc{scaled}.compare.labels.txt"),
    threads: 1
    resources:
        mem_mb=3000,
        time=30,
        partition="high2",
    conda: "conf/env/sourmash.yml"
    log:  os.path.join(logs_dir, "compare", "{basename}.{moltype}.k{k}-sc{scaled}.compare.log")
    benchmark:  os.path.join(logs_dir, "compare", "{basename}.{moltype}.k{k}-sc{scaled}.compare.benchmark")
    shell:
        """
        sourmash compare {input.zip} --ksize {wildcards.k} --avg-containment --ani -o {output.np} --csv {output.csv} 2> {log}
        """

rule plot_aai:
    input: 
        np=os.path.join(out_dir, "{basename}.{moltype}.k{k}-sc{scaled}.compare"),
        labels=os.path.join(out_dir, "{basename}.{moltype}.k{k}-sc{scaled}.compare.labels.txt"),
    output:
        png=os.path.join(out_dir, "plots", "{basename}.{moltype}.k{k}-sc{scaled}.compare.matrix.png"),
    params:
        plot_dir = os.path.join(out_dir, "plots"),
    threads: 1
    resources:
        mem_mb=3000,
        time=30,
        partition="high2",
    conda: "conf/env/sourmash.yml"
    log:  os.path.join(logs_dir, "plot", "{basename}.{moltype}.k{k}-sc{scaled}.plot.log")
    benchmark:  os.path.join(logs_dir, "plot", "{basename}.{moltype}.k{k}-sc{scaled}.plot.benchmark")
    shell:
        """
        mkdir -p {params.plot_dir}
        sourmash plot {input.np} --output-dir {params.plot_dir} --labels 2> {log}
        """

rule plot_pocp_aai_comparison:
    input: 
        compare_csv=os.path.join(out_dir, "{basename}.{moltype}.k{k}-sc{scaled}.compare.csv"),
        pocp_table=POCP_TABLE,
    output:
        comparison_csv=os.path.join(out_dir, "pocp-compare", "{basename}.{moltype}.k{k}-sc{scaled}.pocp-compare.csv"),
        comparison_plot=os.path.join(out_dir, "pocp-compare", "{basename}.{moltype}.k{k}-sc{scaled}.pocp-compare.png"),
    threads: 1
    resources:
        mem_mb=3000,
        time=30,
        partition="high2",
    conda: "conf/env/sourmash.yml"
    log:  os.path.join(logs_dir, "pocp-compare", "{basename}.{moltype}.k{k}-sc{scaled}.pocp-compare.log")
    benchmark:  os.path.join(logs_dir, "pocp-compare", "{basename}.{moltype}.k{k}-sc{scaled}.pocp-compare.benchmark")
    shell:
        """
        python compare-aai-pocp.py --sourmash-compare-csv {input.compare_csv} \
                                   --pocp-table {input.pocp_table} \
                                   --output-comparison-csv {output.comparison_csv} \
                                   --output-plot {output.comparison_plot} 2> {log}
        """