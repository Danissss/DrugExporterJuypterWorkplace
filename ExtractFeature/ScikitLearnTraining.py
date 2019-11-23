# Some guide
# Random Forest
# Feature selection
# Evaluation
import sys, os, csv
import pandas as pd 
import numpy as np
from sklearn.model_selection import cross_val_score
from sklearn.datasets import make_blobs
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.tree import DecisionTreeClassifier
# Removing features with low variance (automatically remove zero and certain threshold)
from sklearn.feature_selection import VarianceThreshold
# Univariate feature selection
from sklearn.datasets import load_iris
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import GenericUnivariateSelect
from sklearn.feature_selection import chi2
# Tree based feature elimination
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.datasets import load_iris
from sklearn.feature_selection import SelectFromModel
# recursive feature elimination
from sklearn.svm import SVC
from sklearn.datasets import load_digits
from sklearn.feature_selection import RFE

from sklearn.ensemble import RandomForestClassifier

from sklearn.model_selection import GridSearchCV # search best hyper-parameters
from sklearn.model_selection import train_test_split

# Load the digits dataset
# digits = load_digits()
# X = digits.images.reshape((len(digits.images), -1))
# y = digits.target

# # Create the RFE object and rank each pixel
# svc = SVC(kernel="linear", C=1)
# rfe = RFE(estimator=svc, n_features_to_select=1, step=1)
# rfe.fit(X, y)
# ranking = rfe.ranking_.reshape(digits.images[0].shape)


# Regulation can reduce the chance of overfitting 
# regularize with L2, early stopping or dropout.


# https://scikit-learn.org/stable/modules/model_evaluation.html



def LoadDataFromCSV(csvFile):
	"""
	Get all descriptor name in list
	Args:
		
	Returns:
		two list: data and class 
	Raise:
		Exceptions
	"""
	data = pd.read_csv(csvFile)
	x = data.iloc[:, 0:315]  # select columns 1 through end
	y = data.iloc[:, 315]   # select column 0, the stock price

	return x,y


def VarianceFeatureSelect(data):
	"""
	Get all descriptor name in list
	Args:
		
	Returns:
		void
	Raise:
		Exceptions
	"""
	sel = VarianceThreshold(threshold=(.8 * (1 - .8)))
	sel.fit_transform(data)


def UnivariateFeatureSelection(data,y,k):
	"""
	Get all descriptor name in list
	Args:
		
	Returns:
		panda data frame: the best value
	Raise:
		Exceptions
	"""

	X_new = SelectKBest(chi2, k).fit_transform(data, y)
	return X_new


def TreeBasedFeatureSelection(data,y):
	"""
	Get all descriptor name in list
	Args:
		param1 (panda frame): attribute value
		param2 (panda frame): class
	Returns:
		panda frame: the transformed data 
	Raise:
		Exceptions
	"""

	# n_estimator = 50
	clf = ExtraTreesClassifier()
	clf = clf.fit(data, y.values.ravel())

	model = SelectFromModel(clf, prefit=True)
	transformed = model.transform(data)
	selected_attribute = data.columns[model.get_support()]
	# print(transformed)
	# transformed_dataframe = pd.DataFrame(data=transformed[1:,1:],index=transformed[1:,0],columns=transformed[0,1:]) 
	# print(transformed_dataframe)

	# print(type(selected_attribute))
	selected_attribute_list = []
	for i in selected_attribute:
		selected_attribute_list.append(i)
	matrix = pd.DataFrame(data=transformed, columns=selected_attribute_list)

	return matrix


def FeatureSelection(data):



	return None


# def CountClass(y):

# 	positive = 0
# 	negative = 0
# 	for i in y:
# 		if i == "none-substrate":



def Train(data,y):


	clf = RandomForestClassifier(n_estimators=100,random_state=0)

	X_train, X_test, y_train, y_test = train_test_split(data, y, test_size=0.2, random_state=0)
	# sys.exit(0)
	# X_train.astype('float64')
	# X_train = X_train.fillna(X_train.median(axis=0))
	clf.fit(X_train, y_train)
	scores = cross_val_score(clf, X_train, y_train, cv=10)
	print("Accuracy: %0.2f (+/- %0.2f)" % (scores.mean(), scores.std() * 2))

	return clf




def main():

	x, y = LoadDataFromCSV("/Users/xuan/Desktop/RDKITWorkPlace/Dataset/Multidrug_resistance_protein_1_non_duplicate_substrate_value.csv")
	
	# CountClass(y)
	# classification without feature selection
	Train(x,y)

	# classification with featuer selection
	new_x = TreeBasedFeatureSelection(x,y)
	Train(new_x,y)

	


if __name__ == '__main__':
	main()



















