configfile: "config.yaml"

# melt each GSM into plus and minus strands in sample table.
import pandas as pd
from numpy import nan

stbl = pd.read_table(config['samples'], dtype = str)
gs = stbl.loc[stbl['sample'].str.startswith('GS')]
gs = gs.assign(plus = pd.Series(['plus'] * gs.shape[0]))
gs = gs.assign(minus = pd.Series(['minus'] * gs.shape[0]))
value_vars = ['plus', 'minus']
other_vars = gs.columns.difference(value_vars)
gs = pd.melt(gs, id_vars = other_vars, value_vars = value_vars, var_name = 'strand', ignore_index = False).sort_index().drop('value', axis = 1)
gs = gs.assign(label_onPlot = pd.Series(gs['label'] + '_' + gs['strand'] + '_' + gs['cellline']))
gs = gs.reindex(sorted(gs.columns), axis = 1)
encode = stbl.loc[stbl['sample'].str.startswith('ENC')]
encode = encode.assign(strand = pd.Series([nan] * encode.shape[0]))
encode = encode.assign(label_onPlot = pd.Series(encode['label'])) 
encode = encode.reindex(sorted(encode.columns), axis = 1)
stbl = pd.concat([gs, encode], axis = 0, ignore_index = True)
stbl = stbl.assign(index = pd.Series((stbl.index + 1).astype(str))) # the unique file name, as str. wildcard captured is a string not number!! 

# define global regex 
from re import escape
wildcard_constraints:
	index = '|'.join([escape(x) for x in stbl['index']]),
	sample = '|'.join([escape(x) for x in stbl['sample']]),
	label = '|'.join([escape(x) for x in stbl['label']]),
	genome = '|'.join([escape(x) for x in set(stbl['genome'])]),
	strand = '|'.join([escape(x) for x in set(x for x in stbl['strand'] if pd.notna(x))])


############# 1.download and renaming (hard-coding GS ENC) 
acceptName = { "minus": '*MINUS*,*minus*', 'plus': '*PLUS*,*plus*'}
rule download_bw:
	output:
		'data/{index}.bigWig'
	params:
		ftp = lambda wc: stbl.loc[stbl['index'] == wc.index, 'ftp'].tolist(),
		accept = lambda wc: stbl.loc[stbl['index'] == wc.index, 'strand'].replace(acceptName).tolist(),
		sample = lambda wc: stbl.loc[stbl['index'] == wc.index,'sample'].tolist()
	shell:
		"""
		if [[ {params.sample} =~ 'GS' ]]
		then 
        	wget --no-remove-listing -r -nd -A '{params.accept}' -O {output} {params.ftp}
        	wget --no-remove-listing -r -nd -A '{params.accept}' -O {output} {params.ftp}
		elif [[ {params.sample} =~ 'ENC' ]]
		then
			wget --no-remove-listing -r -nd -A '{params.sample}' -O {output} {params.ftp}
		else
			echo 'illigal sample pattern {params.sample}'
		fi
		"""
        

rule download_genomeInfo:
	output:
		chainFile = 'data/hg19tohg38',
		chromSize = 'data/hg38.size',
		blacklist = 'data/blacklist_hg38'
	params:
		http_chainFile = config["chainFile"],
		http_chromSize = config["chromSize"],
		http_blacklist = config["blacklistEncode"],
		zip_chainFile =  lambda wildcards, output: output.chainFile +'.gz',
		zip_blacklist =  lambda wildcards, output: output.blacklist +'.gz'
	shell:
		"""
		wget -O {params.zip_chainFile} {params.http_chainFile}
		gunzip -c {params.zip_chainFile} > {output.chainFile}
		rm -f {params.zip_chainFile}

		wget -O {output.chromSize} {params.http_chromSize}
		
		wget -O {params.zip_blacklist} {params.http_blacklist}
		gunzip -c {params.zip_blacklist} > {output.blacklist}
		rm -f {params.zip_blacklist}
		"""


############# 2.liftover of geo data that has been aligned against hg19 
rule hg19tohg38:
	input:
		chainFile = rules.download_genomeInfo.output.chainFile,
		chromSize = rules.download_genomeInfo.output.chromSize,
		hg19_bw = 'data/{index}.bigWig'
	output:
		hg38_bw = 'data/{index}_hg38.bigWig'
	conda:
		'envs/bw_convert.yaml'
	params:
		unmapped = 'data/unmpped',
		hg19_bedgraph = lambda wc: 'data/' + wc.index + '_hg19.bedGraph',
		hg38_bedgraph = lambda wc: 'data/'+ wc.index + '_hg38.bedGraph',
		hg38_sorted_bedgraph = lambda wc: 'data/' + wc.index + '_hg38.sorted.bedGraph',
		hg38_awk_bedgraph = lambda wc: 'data/' + wc.index + '_hg38.awk.bedGraph',
		hg38_fixedInterval_bedgraph = lambda wc: 'data/' + wc.index + '_hg38.fixedInterval.bedGraph', 
		genome = lambda wc: stbl.loc[stbl['index'] == wc.index, 'genome'].tolist() 
	shell:
		"""
		if [ {params.genome} == 'hg19' ]
		then
			# https://www.biostars.org/p/81185/#476941
			bigWigToBedGraph {input.hg19_bw} {params.hg19_bedgraph}
			liftOver {params.hg19_bedgraph} {input.chainFile} {params.hg38_bedgraph} {params.unmapped}
			sort -k1,1 -k2,2n {params.hg38_bedgraph} > {params.hg38_sorted_bedgraph}
			awk -vOFS="\\t" "{{ print \$1, \$2, \$3, \\".\\", \$4 }}" {params.hg38_sorted_bedgraph} > {params.hg38_awk_bedgraph}
			bedops --partition {params.hg38_awk_bedgraph} | bedmap --echo --mean --delim "\\t" - {params.hg38_awk_bedgraph} >  {params.hg38_fixedInterval_bedgraph}
			bedGraphToBigWig {params.hg38_fixedInterval_bedgraph} {input.chromSize} {output.hg38_bw}
			rm -f {input.hg19_bw} {params.hg19_bedgraph} {params.hg38_bedgraph} {params.hg38_sorted_bedgraph} {params.hg38_awk_bedgraph} {params.hg38_fixedInterval_bedgraph} {params.unmapped} 
		elif [ {params.genome} == 'hg38' ]
		then
			mv {input.hg19_bw} {output.hg38_bw}
		else
			echo "genome is {params.genome}"
		fi
		"""

############## 3. deeptools visualization 
rule multiBigWigSummary:
	input: 
		blacklist_bed = rules.download_genomeInfo.output.blacklist,
		bw = expand('data/{index}_hg38.bigWig', index = stbl['index'].tolist())
	output:
		'results/corrMatrix.npz'
	conda:
		'envs/deeptools.yaml'
	params:
		blacklist = 'data/blacklist_hg38',
		labels = ' '.join(stbl['label_onPlot'].tolist()), # space separated
		bs = 10000,
		distBwtBin = 0
	shell:
		"""
		multiBigwigSummary bins -b {input.bw} -o {output} --labels {params.labels} -bs {params.bs} -n {params.distBwtBin} -bl {params.blacklist} -p max
		"""


rule plotCorrelation:
	input:
		rules.multiBigWigSummary.output
	output:
		'results/plot.png'
	conda:
		'envs/deeptools.yaml'
	params:
		c = 'pearson',
		p = 'heatmap',
		colorMap = 'coolwarm'
	shell:
		"""
		plotCorrelation -in {input} --corMethod {params.c} --whatToPlot {params.p} -o {output} --colorMap {params.colorMap} --plotNumbers
		"""
		
############# future implementation: replicate merge
# bigwigMerge (adding score) might not be robust. 
# see https://www.biostars.org/p/209443/
# TODO study how ENCODE merge


############# 
rule all:
	input:
		expand(rules.hg19tohg38.output, index = stbl['index'].tolist()),
		rules.multiBigWigSummary.output,
		rules.plotCorrelation.output
