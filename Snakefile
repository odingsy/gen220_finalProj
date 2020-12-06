configfile: "config.yaml"

# use pandas to load sample table
import pandas as pd
stbl = pd.read_table(config["samples"], dtype = str).set_index("sample", drop=False)
geo = stbl.loc[stbl.index.str.startswith('GS')]
encode = stbl.loc[stbl.index.str.startswith('ENC')]


############# 1.download 
rule download:
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


############# 
rule all:
    input: 
        expand(rules.download.output, sample = geo['sample'])
