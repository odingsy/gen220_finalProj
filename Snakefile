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
        if [ ! -s {output.plus} ]; then
            rm -f {output.plus}
            wget -r -nd -A '*PLUS*,*plus*' -O {output.plus} {params.ftp} # && [[ -s {output.plus} ]]
            sleep 2
        fi

        
        if [ ! -s {output.minus} ]; then
            rm -f {output.minus}
            wget -r -nd -A '*MINUS*,*minus*' -O {output.minus} {params.ftp} # && [[ -s {output.minus} ]]
            sleep 2
        fi
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
    import os
    plus = glob_wildcards('data/{sample}_plus.bigWig').sample
    notplus = list(set(list(geo.index)) - set(plus))
    notplus_plus = [sample + '_plus' for sample in notplus] 
    cmd = 'rm -f ' + ' '.join(expand('data/{yet_to_download}.bigWig', yet_to_download = notplus_plus))
    os.system(cmd)

    minus = glob_wildcards('data/{sample}_minus.bigWig').sample
    notminus = list(set(list(geo.index)) - set(minus))
    notminus_minus = [sample + '_minus' for sample in notminus]
    cmd = 'rm -f ' + ' '.join(expand('data/{yet_to_download}.bigWig', yet_to_download = notminus_minus))
    os.system(cmd)

    return expand('data/{yet_to_download}_{ps}.bigWig', yet_to_download = notplus + notminus, ps = ['plus', 'minus'])



############# 
rule all:
    input:
        BounceBackReDownload
        #expand(rules.download.output, sample = geo['sample'])
