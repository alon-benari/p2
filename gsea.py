import pandas as pd
from itertools import permutations
import numpy as np
import matplotlib.pyplot as plt
import random
import sys


class GSEA():

  def __init__(self,expfile,sampfile,keggfile):
    """
    Initialize the GSEA class with file names

    expfile -  gene expression file
    smplfile -  samples for healthy and sick individuals
    keggfile -  a file with gene set assignments.
    """

    self.expression = pd.read_table(expfile,sep = '\t')
   
    self.samples = pd.read_table(sampfile,sep = '\t',names = ['samples','patient'])
    #
    self.original_gene_set_es = {}
    self.kegg_es ={} # hold the distribution of ES for each gene set
    self.es_dict = {} # a dictionary to hold the ESs 
    self.ratio_dict ={} # a dictionary to hold the ranked ratio scores key by number
    self.kegg_dict = {} # a dictionary  keyed by gene set name and value of gene set
    self.shuffle_dict ={} # a dictionary to hold the shuffled labels
    f = open(keggfile)
    kegg_lines = f.readlines()
    #
    for l in kegg_lines:
        rec = l.replace('\t',' ').strip('\n').split(' ')
        self.kegg_dict[rec[0]] = rec[2:] 
    #
    self.patient = self.samples[self.samples['patient'] == 1].samples.tolist()
    self.healthy = self.samples[self.samples['patient'] == 0].samples.tolist()

    #
    self.set_data() # a function to wrangle the data into a better shape.
    self.original_lbls = self.expression.columns.tolist()  # preserve real labels
    self.shuffle_labels = (lbls for lbls in permutations(self.samples.index)) # use a generator to get an exhaustive set of permutations of patient samples
    #
  def random_walk(self):
    """
    Take a gene set and look at the ranked list that we have and compute if it is toward the  up regulated or under regulated
    """
    #
    for kegg_k, gene_set in [(k,v) for k,v in self.kegg_dict.items()]:  # pick a bunch of genes
     
      reward = np.sqrt((self.expression.shape[0] - len(gene_set))/float(len(gene_set))) # 
      penalty = -np.sqrt(len(gene_set)/float(self.expression.shape[0]-len(gene_set)))
      #      #
      #print(kegg_k,reward,penalty) # just for debug
    
      for key in self.ratio_dict.keys(): # iterate over the shuffle samples set 
        # print('key: ',key)
        trace = []
        for r in self.ratio_dict[key].index: # this one hold 100 ranked set of the shuffle genes and the original at [0]. ratio_dict is indexed by the genes hence the .index method
          if r in gene_set: # check for existance the gene set
            trace.extend([reward])  # get reward if there
          else:
            trace.extend([penalty]) #get penalty if not
        self.es_dict[key] = np.array(trace).cumsum()# traces of Brownian Bridges key by shuffle id 
      self.kegg_es[kegg_k] = self.es_dict
      self.es_dict = {}

  def show_trace(self,gene_name):
    """
    A method to return the traces for the gene set
    """
    self.kegg_es


  def output_es(self):
    """
    A method to output the ES for the gene set before shuffling the sample expression data
    """
    rows = []
    for k in self.kegg_es.keys():
      ES_true = [format(v.max(), '.2f') for k,v in self.kegg_es[k].items()][0]
      rows.append([k,ES_true])
    #
    output = pd.DataFrame(rows, columns = ['gene_set','es'])
    output.set_index('gene_set',inplace = True)
    out2file = output['es'].sort_values(ascending=False)
    out2file.to_csv('kegg_enrichment_scores.txt',sep = '\t')
    return out2file

  def p_value(self,cutoff=0.05):
    """
    A method to return the gene sets where enrichment 
    """
    p_val ={}
    for k,v in self.kegg_es.items():
      es_dist = [val.max() for key, val in v.items()]
      p_val[k] = sum(es_dist[1:]>es_dist[0])/float(len(es_dist))
    return  p_val #np.array([v for k,v in p_val.items()])
      

  def true(self):
    """
    A method to return the ranking of the original matrix
    """
    return self.ratio_dict[0].index


  def true_rank(self):
    """
    A method to return  the true ranked ratio from a data frame  before shuffling  the columns
    input -  a dataframe of samples and genes expressions with the CORRECT REAL labeling
    output -  dictionary of ranked ratios
    """
    set1 = self.expression[self.patient]
    set2 = self.expression[self.healthy]
    # ratio = set1.mean(axis=1)/set2.mean(axis=1) # compute ratio
    ratio = set1.mean(axis=1)-set2.mean(axis=1)
    self.ratio_dict[0] = ratio.sort_values(ascending = False) # store the ranked result in key 0 for ES_true


  def compute_ranks(self):
    """
    A method to return ranked means ratio  of the said sample
    output - 
      returns a dictionary of 1000 ranked ratios
    """
    
    # relabel the sample:
    for k in range(1,100):
      temp_df = self.expression # this is the expression set
      #l = next(self.shuffle_labels) # pull a fresh set of shuffled labels
      l = random.sample(list(self.samples.index),len(list(self.samples.index)))
      #
      # check that we dont get the original column sequence
      if (''.join(list(l)) != ''.join(self.original_lbls)):
        self.shuffle_dict[k] = l    # keep track of columns
        temp_df.columns = l  # replace shuffle label over data set, this is where I assign new columns from l
        #       
        ratio = temp_df[self.patient].mean(axis=1)-temp_df[self.healthy].mean(axis=1) # compute log ratio 
        self.ratio_dict[k] = ratio.sort_values(ascending = False) # store the ranked result
      else:
        print('woops no can do, move to next  combination, this is the correct one.')

     
  
  def set_data(self):
    """
    A method to set the data in the right  format
    
    It basically returns two data frames with the samples as indices
    no input or output.
    """
    # take care of samples
    patients = self.samples.iloc[:,1].tolist()
    samples = self.samples.iloc[:,0].tolist()
    self.samples  = pd.DataFrame(patients,index = samples,columns = ['patient']) # indexed by sample
    #
    self.expression = self.expression.set_index('SYMBOL')
    # # take care of expression data
    # self.data = pd.merge(self.expression,self.samples,left_index = True,right_index=True) # merged data sets
    # #pd.merge(df1,df2,how = 'left',left_index=True,right_index=True) # do a left join

  def kegg_genes(self):
    """
    A method to return the number of unique genes in the  KEGG set
    """
    gene_pool = []
    uniq = {} # a dictionary to hold the counts of genes.
    for k,v in self.kegg_dict.items():
      gene_pool.extend(v)
    for u in set(gene_pool):
      uniq[u] = len(np.where(np.array(gene_pool) == u)[0])
    
    kegg_df = pd.DataFrame.from_dict(uniq,orient = 'index')
    kegg_df.columns = ['counts']
    return kegg_df

  def  max_es(self):
    """
    A method to return the gene with the maximal enrichment set
    """
    es_original = {k:v[0].max() for k,v in self.kegg_es.items()}
    es = pd.DataFrame.from_dict(es_original,orient = 'index')
    es.columns = ['ES']

    return es.sort_values('ES')

 
    
#### run 
def main():
  expfile = sys.argv[1] #'GSE25628_filtered_expression.txt'
  sampfile =  sys.argv[2] #'GSE25628_samples.txt'
  keggfile = sys.argv[3] #'c2.cp.kegg.v6.2.symbols.filtered.gmt'

  gsea = GSEA(expfile,sampfile,keggfile)
  gsea.true_rank()
  gsea.compute_ranks()  
  gsea.random_walk()
  o = gsea.output_es() # output toa file
  pval = gsea.p_value()
  print('Number of Significant pathways=',np.array(np.array([v for k,v in pval.items()])< 0.05/145.0).sum())
  max_es = gsea.max_es()


if __name__ == '__main__':
  main()
  # python3 gsea.py  GSE25628_filtered_expression.txt  GSE25628_samples.txt c2.cp.kegg.v6.2.symbols.filtered.gmt
