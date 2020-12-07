configfile: "config.yaml"

# use pandas to load sample table
import pandas as pd
stbl = pd.read_table(config["samples"], dtype = str).set_index("sample", drop=False)
geo = stbl.loc[stbl.index.str.startswith('GS')]
encode = stbl.loc[stbl.index.str.startswith('ENC')]


############# 1.download 
checkpoint download:
    output:
        plus='data/{sample}_plus.bigWig',
        minus='data/{sample}_minus.bigWig'
    params:
        ftp=lambda wc: geo.loc[geo['sample'] == wc.sample, 'ftp'].item(),
        # retryConn='--retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 3 --no-dns-cache'
        # https://superuser.com/questions/493640/how-to-retry-connections-with-wget
    shell:
        """
        sleep 2
        wget -r -nd -A '*PLUS*,*plus*' -O {output.plus} {params.ftp} # && [[ -s {output.plus} ]]
        sleep 2
        wget -r -nd -A '*MINUS*,*minus*' -O {output.minus} {params.ftp} # && [[ -s {output.minus} ]]
        sleep 2
        """
        

############# and renaming (GS: label_plus.bw, ENC:label.bw) 


############# 2.liftover 
# https://www.biostars.org/p/81185/#476941

############# future implementation: replicate merge
# bigwigMerge (adding score) might not be robust. 
# see https://www.biostars.org/p/209443/
# TODO study how ENCODE merge


#############  

def BounceBackReDownload(wildcards):
    import os, re
    from itertools import product

    sampdir = 'data'
    existingSamp = []

    # get rid of empty file. 
    for samp in os.scandir(sampdir):
        if os.path.getsize(samp) == 0:
            os.remove(samp)
        else:
            existingSamp.append(samp.path)
    allSamp = [samp + pm + '.bigWig' for samp in geo.index for pm in ['_plus', '_minus']]
    allSamp = [os.path.join(sampdir, samp) for samp in allSamp]
    redownloadSamp = list(set(allSamp) - set(existingSamp))
    # samples = [ re.search('/(.+?)_', samp)[1] for samp in redownloadSamp ]
    # pms = [ re.search('_(.+?)\.', samp)[1] for samp in redownloadSamp ]
    return redownloadSamp



############# 
rule all:
    input:
        BounceBackReDownload
        #expand(rules.download.output, sample = geo['sample'])
