__author__ = 'lucas Ferreira'
import simplehmm
import re
treinamento=[]
arquivo_motivos=open('teste.mtv','r').read()
for motivo in arquivo_motivos.split('\n'):
    motivo_treino=[]
    contador=0
    for x in re.findall('\w',motivo):
        motivo_treino.append([str(contador),x])
        contador+=1
    treinamento.append(motivo_treino)

#print treinamento
cirRNA_hmm=simplehmm.hmm('circ_rna primeiros testes',['0','1','2','3','4','5','6','7','8','9','10'], ['A','C','T','G'])
cirRNA_hmm.train(treinamento, smoothing='absdiscount')
cirRNA_hmm.print_hmm()