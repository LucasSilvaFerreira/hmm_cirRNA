__author__ = 'lucas Ferreira'
import  simplehmm
import re

cirRNA_hmm=simplehmm.hmm('circ_rna primeiros testes',['dummy'], ['dummy'])
cirRNA_hmm.load_hmm('circRNA_FRAGMENTADA.hmm')


arquivo = open('original_sem_deletar_os_com_gap_introns_ciRNAs.txt','r').read()
for seq in arquivo.split('>'):
    if len(seq)>0:
        header= seq.split('\n')[0]
        position= header.split('_')[6].split('-')
        strand= header.split('_')[7]
        branch= int(re.sub(';','',header.split('_')[8]))
        print branch,header
        #print strand,position
    for seq_bases in re.split('hg19.*\n' ,seq):
        sequencia_completa = re.sub('\r|\n','',seq_bases)
        #print sequencia_completa
    contador=0

    probabilidades=[]
    for x in range(0,len(sequencia_completa)):
        query=[]
        if x < len(sequencia_completa)-11:
            for y in re.findall('\w',sequencia_completa[x:x+11]):
                query.append(y)
            probabilidades.append(cirRNA_hmm.viterbi(query)[1])

    try:
        #print len (probabilidades)
        print sorted(probabilidades)[len(probabilidades)-1]
        print probabilidades.index(sorted(probabilidades)[len(probabilidades)-1])
        print int(position[0])+int(probabilidades.index(sorted(probabilidades)[len(probabilidades)-1])),branch
    except:
        print 'zero'
