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
        sequencia_completa=sequencia_completa[::-1]
    contador=0
    sequencia_motivo=[]
    print (sequencia_completa)
    probabilidades=[]
    for x in range(0,len(sequencia_completa)):
        query=[]
        if x < len(sequencia_completa)-11:
            for y in re.findall('\w',sequencia_completa[x:x+11]):
                query.append(y)


            if query[0]=='A'or query[0]=='C':
                if query[7]=='C' or query[7]=='G':
                    if query[9]=='C' or query[9]=='T':
                        probabilidades.append([cirRNA_hmm.viterbi(query)[1],sequencia_completa[x:x+11],x])
                        sequencia_motivo.append(sequencia_completa[x:x+11])


            else:  # necessario para nao danificar o contador de posicao dentro do vetor
                probabilidades.append(-1)
                sequencia_motivo.append(-1)
    try:
        #print len (probabilidades)
        #print sorted(probabilidades)[len(probabilidades)-1]
        for i in range(1,10):
            print'-----------------------------------'
            print sorted(probabilidades)[len(probabilidades)-i]
            #print probabilidades.index(sorted(probabilidades)[len(probabilidades)-i])
            #print sequencia_motivo[probabilidades.index(sorted(probabilidades)[len(probabilidades)-i])]
        print int(position[0])+int(probabilidades.index(sorted(probabilidades)[len(probabilidades)-1])),branch
    except:
        print 'zero'
