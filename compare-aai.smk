import pandas as pd

out_dir = "output.pocp"
logs_dir = os.path.join(out_dir, 'logs')
basename = 'brady'

rule all:
    input: 
        expand(os.path.join(out_dir, f"{basename}.{{moltype}}.sc{{scaled}}.zip"), moltype=["protein"], scaled=[5,10]),
        expand(os.path.join(out_dir, f"{basename}.{{moltype}}.sc{{scaled}}.compare.csv"), moltype=["protein"], scaled=[5,10]),
        expand(os.path.join(out_dir, "plots", f"{basename}.{{moltype}}.sc{{scaled}}.compare.matrix.png"), moltype=["protein"], scaled=[5,10])

paramD = {"protein": "protein,k=10,scaled=1,abund"}
rule sketch_fromfile:
    input: 
        fromfile=os.path.join(out_dir, "{basename}.fromfile.csv"),
    output: os.path.join(out_dir, "{basename}.{moltype}.zip")
    params:
        lambda w: paramD[w.moltype]
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
        sourmash sketch fromfile {input.fromfile} -p {params} -o {output} 2> {log}
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
        csv=os.path.join(out_dir, "{basename}.{moltype}.sc{scaled}.compare.csv"),
        np=os.path.join(out_dir, "{basename}.{moltype}.sc{scaled}.compare"),
        labels=os.path.join(out_dir, "{basename}.{moltype}.sc{scaled}.compare.labels.txt"),
    params:
        lambda w: paramD[w.moltype]
    threads: 1
    resources:
        mem_mb=3000,
        time=30,
        partition="high2",
    conda: "conf/env/sourmash.yml"
    log:  os.path.join(logs_dir, "compare", "{basename}.{moltype}.sc{scaled}.compare.log")
    benchmark:  os.path.join(logs_dir, "compare", "{basename}.{moltype}.sc{scaled}.compare.benchmark")
    shell:
        """
        sourmash compare {input.zip} --avg-containment --ani -o {output.np} --csv {output.csv} 2> {log}
        """

rule plot_aai:
    input: 
        np=os.path.join(out_dir, "{basename}.{moltype}.sc{scaled}.compare"),
        labels=os.path.join(out_dir, "{basename}.{moltype}.sc{scaled}.compare.labels.txt"),
    output:
        png=os.path.join(out_dir, "plots", "{basename}.{moltype}.sc{scaled}.compare.matrix.png"),
    params:
        plot_dir = os.path.join(out_dir, "plots"),
    threads: 1
    resources:
        mem_mb=3000,
        time=30,
        partition="high2",
    conda: "conf/env/sourmash.yml"
    log:  os.path.join(logs_dir, "plot", "{basename}.{moltype}.sc{scaled}.plot.log")
    benchmark:  os.path.join(logs_dir, "plot", "{basename}.{moltype}.sc{scaled}.plot.benchmark")
    shell:
        """
        mkdir -p {params.plot_dir}
        sourmash plot {input.np} --output-dir {params.plot_dir} --labels 2> {log}
        """