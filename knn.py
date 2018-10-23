import numpy as np
import pandas as  pd
import sys


class  KNN():
  """
  A class to implement KNN and report performance of classifier.
  """
  def __init__(self,k,p,expfile,sampfile):
    """
    A method to initialize a KNN with the number of k nearest number
     k - k nearest nighbor
     p - fraction of positive patient needed to accept 'positive'
     expfile -  exrpession file
     sampfile -  file of patients and diagnosis
    """
  
    self.k = k
    self.p = p
    self.expression =pd.read_table(expfile,sep = '\t')
    self.samples = pd.read_table(sampfile,sep = '\t',names = ['samples','patient'])
    self.set_data() # a functionto wrangle the data a bit.

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
    # take care of expression data
    cols = self.expression.SYMBOL.tolist()  # set new column names to transposed expression_data 
   
    new_exp = self.expression.T.ix[1:,:] # transpose
    new_exp.columns = cols
    self.expression = new_exp                    # add columns
    self.data = pd.merge(self.expression,self.samples,left_index = True,right_index=True) # merged data sets
    #pd.merge(df1,df2,how = 'left',left_index=True,right_index=True) # do a left join
  
  def loo(self):
    """
    A method to retunr a training set and a test set.
    get the indices for leave one out for all sets
    return-
      a list of lists [[train_index], [test_index]]
    """
    loo = list()
    for i in range(self.data.shape[0]):
      train_index = [i for i in range(self.data.shape[0])]
      train_index.pop(i)
      test_index = i
      loo.append([train_index,test_index])

    return (loo)

  def train(self,X,y):
    """
    A method to initialize the model
    parameters-

    X ,y  -  training set and training labels repsectively
    """
    self.X_train = X
    self.y_train = y
    self.class_labels = np.unique(self.y_train)

   
  
  def predict(self,X,y):
    """
    A method to predict the class label of a sample
    parameters - 
    X,y - input a test set and test label as numpy arrays or pandas data frame/ series
    output - 
      dictionary :
        sample id,y_pred
    """
    self.X_test = X
    self.y_test = y
    d = []
    for i in range(self.X_train.shape[0]):
      d.append(self.get_distance(self.X_train.ix[i,:])) # hold all distances
    sorted = np.argsort(d)
    k_indices = np.argsort(d)[:self.k] # get indices with lowest distances
    predictions = self.y_train[k_indices]
    unique, counts =  np.unique(predictions,return_counts=True)

    if (np.where(predictions ==1)[0].shape[0]) >self.p*self.k:
      y_pred = 1
    else:
      y_pred=0
      # {'sample':X_test.name,'d':d,'k_ix':k_indices,'pred':predictions,
      #       'counts':counts,'uniq':unique,'y_pred':y_pred,
      #       'y_test':self.y_test,'y_train':self.y_train,
      #       'sorted':sorted}
    return {'sample':self.X_test.name,
            'y_pred':y_pred, 
            'y_test':self.y_test}



  def get_distance(self,row_vector):
    """
    A method to compute the distance between  X_test and X_train and return the euclidean distance
    input - 
      a pandas Series or numpy array
    output - 
      euclidean distance
    """
    d = row_vector-self.X_test
    
    return np.sqrt(np.dot(d,d)) # return the euclidean distance
    
  def cm(self,cm):
    """
     A method to returb a confusion matrix, TP rate and FP rate
     input - 
        a data frame y_pred, and y_true  from predict method
      output -
        dictionary:
          cm: confusion matrix
          sens: senstivity ;tp/(tp+fn)
          spec: specificity; tn/(tn+fp)
    """
    tp = cm[(cm.y_true ==1) & ( cm.y_pred ==1 )].shape[0]
    fn = cm[(cm.y_true == 1) & (cm.y_pred == 0)].shape[0]
    fp = cm[(cm.y_true ==0) & (cm.y_pred == 1)].shape[0]
    tn = cm[(cm.y_true == 0) &( cm.y_pred == 0)].shape[0]
    sens = format(tp/float(tp+fn),'.2f')
    spec = format(tn/float(tn+fp),'.2f')
    return sens, spec

##### run ######
def main():
  expfile =  sys.argv[1] #'GSE25628_filtered_expression.txt' #
  sampfile = sys.argv[2] #'GSE25628_samples.txt' #
  k = int(sys.argv[3])#5 # 
  p = float(sys.argv[4]) # 0.5 
#   #
  knn = KNN(k,p,expfile,sampfile)  # instantiate the object
  y_pred = list()
  y_true = list()
  sample = list()
  for train_index, test_index in knn.loo():
    X_train = knn.data.ix[train_index,:-1]
    y_train = knn.data.ix[train_index].patient
    X_test = knn.data.ix[test_index,:-1]
    y_test = knn.data.ix[test_index].patient
        
    knn.train(X_train,y_train)
    knn_predictions = knn.predict(X_test,y_test)
    y_pred.append(knn_predictions['y_pred'])
    sample.append(knn_predictions['sample'])
    y_true.append(y_test)

  cm = (pd.DataFrame({'sample':sample,'y_pred':y_pred,'y_true':y_true}))
  #
  # write to file
  cm[['sample','y_pred']].to_csv('sample_assignments.txt',sep = '\t',header=False)
  # output to crt
  sens, spec = knn.cm(cm)
  print('sensitivity: ',sens)
  print('specificity: ',spec)

if __name__ =='__main__':
    main()