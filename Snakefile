configfile: "config.yaml"

# use pandas to load sample table
import pandas as pd
stbl = pd.read_table(config["samples"], dtype = str).set_index("sample", drop=False)

import re
wildcard_constraints:
	sample='|'.join([re.escape(x) for x in stbl['sample']]),
	label='|'.join([re.escape(x) for x in stbl['label']]),
	genome='|'.join([re.escape(x) for x in set(stbl['genome'])]),
	pm='|'.join(re.escape(x) for x in ['plus', 'minus'])


############# 1.download and renaming (GS: GSMxxx_plus.bw, ENC:label.bw) 
rule download_geo:
    output:
        plus='data/{sample}_plus.bigWig',
        minus='data/{sample}_minus.bigWig'
    params:
        ftp=lambda wc: stbl.loc[stbl['sample'] == wc.sample, 'ftp'].item()
    shell:
        """
        wget --no-remove-listing -r -nd -A '*PLUS*,*plus*' -O {output.plus} {params.ftp}
        wget --no-remove-listing -r -nd -A '*MINUS*,*minus*' -O {output.minus} {params.ftp}
        """
        
rule download_encode:
	output:
		'data/{label}.bigWig'
	params:
		ftp=lambda wc:stbl.loc[stbl['label']==wc.label,'ftp'].item(),
		sample=lambda wc:stbl.loc[stbl['label']==wc.label,'sample'].item()
	shell:
		"""
		wget --no-remove-listing -r -nd -A '{params.sample}' -O {output} {params.ftp}
		"""

rule download_genomeInfo:
	output:
		chainFile='data/hg19tohg38',
		chromSize='data/hg38.size'
	params:
		http_chainFile = config["chainFile"],
		http_chromSize = config["chromSize"],
		unzip_output =  lambda wildcards, output: output.chainFile +'.gz'
	shell:
		"""
		wget -O {params.unzip_output} {params.http_chainFile}
		gunzip -c {params.unzip_output} > {output.chainFile}
		wget -O {output.chromSize} {params.http_chromSize}
		"""
############# 2.liftover on geo data aligned against hg19 
# https://www.biostars.org/p/81185/#476941
rule hg19tohg38:
	input:
		hg19_bw = 'data/{sample}_{pm}.bigWig'
	output:
		hg19_bedgraph = temp('data/{sample}_{pm}_hg19.bedGraph'),
		hg38_bedgraph = temp('data/{sample}_{pm}_hg38.bedGraph'),
		hg38_sorted_bedgraph = temp('data/{sample}_{pm}_hg38.sorted.bedGraph'),
		hg38_awk_bedgraph = temp('data/{sample}_{pm}_hg38.awk.bedGraph'),
		hg38_fixedInterval_bedgraph = temp('data/{sample}_{pm}_hg38.fixedInterval.bedGraph'),
		hg38_bw = 'data/{sample}_{pm}_hg38.bigWig'
	conda:
		'envs/bw_convert.yaml'
	params:
		unmapped='data/unmpped'
	shell:
		"""
		bigWigToBedGraph {input.hg19_bw} {output.hg19_bedgraph}
		liftOver {output.hg19_bedgraph} data/hg19tohg38 {output.hg38_bedgraph} {params.unmapped}
		sort -k1,1 -k2,2n {output.hg38_bedgraph} > {output.hg38_sorted_bedgraph}
		awk -vOFS="\\t" "{{ print \$1, \$2, \$3, \\".\\", \$4 }}" {output.hg38_sorted_bedgraph} > {output.hg38_awk_bedgraph}
		bedops --partition {output.hg38_awk_bedgraph} | bedmap --echo --mean --delim "\\t" - {output.hg38_awk_bedgraph} >  {output.hg38_fixedInterval_bedgraph}
		bedGraphToBigWig {output.hg38_fixedInterval_bedgraph} data/hg38.size {output.hg38_bw}
		"""	



############# future implementation: replicate merge
# bigwigMerge (adding score) might not be robust. 
# see https://www.biostars.org/p/209443/
# TODO study how ENCODE merge


############# 
rule all:
	input:
		expand(rules.download_geo.output, sample = stbl.loc[stbl.index.str.startswith('GS')]['sample']),
		expand(rules.download_encode.output, label = stbl.loc[stbl.index.str.startswith('ENC')]['label']),
		rules.download_genomeInfo.output, 
		expand(rules.hg19tohg38.output, sample = stbl.loc[stbl['genome'] == 'hg19', 'sample'], pm = ['plus', 'minus'])
